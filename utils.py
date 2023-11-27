import requests
import pandas as pd
from vcf_parser import VCFParser


def parse_vcf_file(file):
    """
    Parse .vcf file into a pd.DataFrame containing defined list of columns. Columns are:
    'ID' SPDI notation
    'CHROM' Chromosome on which variant is located
    'POS' Position of the variant
    'REF' Reference allele
    'ALT' Alternative allele
    'TC' Total coverage at this locus
    'TR' Total number of reads containing this variant
    'AR' Total number of reads containing alternative variant
    'WS' Starting position of calling window
    'WE' End position of calling window
    :param file: .vcf file
    :return: pd.DataFrame with information from .vcf file
    """
    my_parser = VCFParser(infile=file, split_variants=True, check_info=True)

    variants = []
    for variant in my_parser:
        variants.append(variant)
        # Total coverage
        variant["TC"] = int(variant["info_dict"]["TC"][0])
        # Reads containing allele
        variant["TR"] = int(variant["info_dict"]["TR"][0])
        # Reads containing alternate allele
        variant["AR"] = variant["TC"] - variant["TR"]
        # Starting position of calling window
        variant["WS"] = variant["info_dict"]["WS"][0]
        # Ending position of calling window
        variant["WE"] = variant["info_dict"]["WE"][0]

        # This is a REST-style regions to query EnsEMBL by
        # variant['ID'] = f"{variant['CHROM']}:{variant['WS']}:{variant['WE']}/{variant['REF']}"
        variant[
            "ID"
        ] = f"{variant['CHROM']} {variant['POS']} . {variant['REF']} {variant['ALT']} . . ."
    variants = pd.DataFrame(variants)
    variants = variants.loc[
        :, ["ID", "CHROM", "POS", "REF", "ALT", "TC", "TR", "AR", "WS", "WE"]
    ]
    return variants


def get_response(request, species, chunk_size=100):
    """
    Get response from EnsEMBL REST API for a certain variant using SPDI notation list as input
    :param request: List of variants in SPDI notation
    :param species: Species for which query is run
    :param chunk_size: Chunk to which list is split in separate blocks to be sent by post request
    :return: List of responses
    """
    url = f"https://rest.ensembl.org/vep/{species}/region/"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    responses = []
    for i in range(0, len(request), chunk_size):
        chunk = request[i : i + chunk_size]
        chunk = "{" + f""""variants": ["{'", "'.join(chunk)}"]""" + "}"
        response = requests.post(url, headers=headers, data=chunk)
        if response:
            response = response.json()
            responses.extend(response)
            continue

    return responses


def process_response(data):
    """
    Process dict response from EnsEMBL REST API into pd.DataFrame with the following columns:
    "ID" - SDPI ID which was in the request
    "consequence_group" - group of the consequences of the variant
    "consequence_terms" - consequence of the variant
    "gene_symbol" - gene symbol
    "allele" - allele of the variant
    :param data: Response from EnsEMBL REST API, dictionary.
    :return: pd.DataFrame with a data from response.
    """
    input = {}
    output = {}
    for response in data:

        input_id = response["input"]

        for consequences in [
            "intergenic_consequences",
            "transcript_consequences",
            "regulatory_feature_consequences",
        ]:
            if consequences in response.keys():
                input[consequences] = response[consequences]

        for key, value in input.items():
            output[input_id, key] = {}
            for d in value:
                output[input_id, key]["allele"] = d["variant_allele"]
                if "gene_symbol" in d.keys():
                    output[input_id, key]["gene_symbol"] = d["gene_symbol"]
                else:
                    output[input_id, key]["gene_symbol"] = ""
                if "consequence_terms" in d.keys():
                    output[input_id, key]["consequence_terms"] = ", ".join(
                        d["consequence_terms"]
                    )
                else:
                    output[input_id, key]["consequence_terms"] = ""

        output = pd.DataFrame.from_dict(output, orient="index")
        output = output.reset_index(names=["ID", "consequences_group"])
    return output

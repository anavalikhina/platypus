import requests
import pandas as pd
from vcf_parser import VCFParser


def parse_vcf_file(file):
    my_parser = VCFParser(infile=file, split_variants=True, check_info=True)

    variants = []
    for variant in my_parser:
        variants.append(variant)
        # Total coverage
        variant['TC'] = int(variant['info_dict']['TC'][0])
        # Reads containing allele
        variant['TR'] = int(variant['info_dict']['TR'][0])
        # Reads containing alternate allele
        variant['AR'] = variant['TC'] - variant['TR']
        # Starting position of calling window
        variant['WS'] = variant['info_dict']['WS'][0]
        # Ending position of calling window
        variant['WE'] = variant['info_dict']['WE'][0]

        # This is a REST-style regions to query EnsEMBL by
        variant['ID'] = f"{variant['CHROM']}:{variant['WS']}:{variant['WE']}/{variant['REF']}"

    variants = pd.DataFrame(variants)
    variants = variants.loc[:, ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'TC', 'TR', 'AR', 'WS', 'WE']]
    return variants


def get_response(request):
    url = f'https://rest.ensembl.org/vep/ornithorhynchus_anatinus/region/'
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    response = requests.post(url, headers=headers, data=request)

    data = response.json()
    return data
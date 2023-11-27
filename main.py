import json
import pandas as pd
from utils import parse_vcf_file, get_response, process_response

with open("parameters.json", "r") as file:
    parameters = json.load(file)

# parsing .vcf file into a dataframe
variants = parse_vcf_file(file=parameters["File"])

# requesting EnsEMBL DB
request = variants["ID"].values.tolist()
response = get_response(
    request, species=parameters["Species"], chunk_size=parameters["Chunk_size"]
)

# processing response into a dataframe
response_dataframe = process_response(response)

# merging processed .vcf with response from EnsEMBL
result = pd.merge(variants, response_dataframe)

# saving file
result.to_csv("result.csv")

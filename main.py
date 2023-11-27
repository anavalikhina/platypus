import json
import pandas as pd
from utils import parse_vcf_file, get_response, process_response

with open("parameters.json", 'r') as file:
    parameters = json.load(file)

variants = parse_vcf_file(file=parameters["File"])

request = variants['ID'].values.tolist()
response = get_response(request, species=parameters["Species"], chunk_size=parameters["Chunk_size"])

response_dataframe = process_response(response)

result = pd.merge(variants, response_dataframe)




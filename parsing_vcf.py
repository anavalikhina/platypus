from utils import parse_vcf_file, get_response


variants = parse_vcf_file(file='test_vcf.vcf')
request = variants['ID'].values.tolist()
# trying piece of request
request = '''{
    "variants": [
        "X:153994586:153994604/G",
        "X:154004992:154005010/G",
        "X:154020074:154020092/C",
        "X:154020104:154020122/C",
        "X:154456737:154456755/A",
        "X:155125424:155125443/AG",
        "X:155127665:155127683/A",
        "X:155233088:155233106/T"]}'''
get_response(request)

#
# if 'transcript_consequences' in data[0].keys():
#     allele_data = data[0]['transcript_consequences'][0]
# else:
#     return '', ''
# if 'gene_symbol' in allele_data.keys():
#     name = allele_data['gene_symbol']
# else:
#     name = ''
# if 'consequence_terms' in allele_data.keys():
#     consequence = allele_data['consequence_terms']
# else:
#     consequence = ''
#



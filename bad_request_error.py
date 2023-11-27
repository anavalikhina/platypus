import requests, sys

# For some reason I have 400 bad request error even with the example code.
# This issue is found only for platypus, but not for human.
server = "https://grch37.rest.ensembl.org"
ext = "/vep/ornithorhynchus_anatinus/region/X:153994586:153994604/G"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
print(repr(decoded))

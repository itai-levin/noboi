import requests
import re
from tqdm import tqdm

METACYC_PROTEINS_PATH = '/Users/Itai/Library/Mobile Documents/com~apple~CloudDocs/from_box/Grad/research/molecule_databases/Metacyc_v26.5/data/protein-seq-ids-unreduced.dat'

uniprot_accessions = []
with open(METACYC_PROTEINS_PATH, 'r') as f:
    for line in f:
        uniprot_accessions.extend(re.findall('(?<=UNIPROT:)[0-9A-Za-z]*', line))


batched_uniprot_accessions = [uniprot_accessions[i:i + 25] for i in range(0, len(uniprot_accessions), 25)]

returned_uniprot = []
for b in tqdm(batched_uniprot_accessions):
    query_string = '+OR+'.join(['accession:{}'.format(s) for s in b])
    res = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={query_string}&format=fasta").text
    returned_uniprot.append(res)


with open('uniprot_chkpoint.txt', 'w') as f:
    for line in returned_uniprot:
        f.write(line)

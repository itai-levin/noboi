import json
import pandas as pd
import requests
import re
from tqdm import tqdm

rhea_2_uniprot = pd.read_csv('rhea/rhea2uniprot_sprot.tsv', sep='\t')

uniprot_accessions = rhea_2_uniprot['ID'].tolist()

#check if its already_queried 
with open('protein_id_to_sequence.json','r') as f:
    uniprot2seq = json.load(f)

uniprot_accessions = [u for u in uniprot_accessions if u not in uniprot2seq.keys()]

print (len(uniprot_accessions), 'uniport accessions to query')

batched_uniprot_accessions = [uniprot_accessions[i:i + 25] for i in range(0, len(uniprot_accessions), 25)]

returned_uniprot = []
for b in tqdm(batched_uniprot_accessions):
    query_string = '+OR+'.join(['accession:{}'.format(s) for s in b])
    res = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={query_string}&format=fasta").text
    returned_uniprot.append(res)


with open('uniprot_chkpoint_for_rhea.txt', 'w') as f:
    for line in returned_uniprot:
        f.write(line)

import json
import pandas as pd
import requests
import re
from tqdm import tqdm

with open('brenda/scraped_brenda_substrate_results.json', 'r') as f:
    parsed_brenda = json.load(f)
for ec in parsed_brenda.keys():
    parsed_brenda[ec] = pd.DataFrame(parsed_brenda[ec])

uniprot_accessions = set()
for ec in parsed_brenda:
    uniprot_accessions=uniprot_accessions.union(set(parsed_brenda[ec]['UNIPROT'].tolist()))
    
ls = [x.replace(' ','').replace(';',',').split(',') for x in list(uniprot_accessions)]
uniprot_accessions = set([x for y in ls for x in y])

#check if its already_queried 
with open('20Nov2023_protein_id_to_sequence.json','r') as f:
    uniprot2seq = json.load(f)

#if os.path.isfile('uniprot_chkpoint_for_brenda.txt'):
#    with open()


uniprot_accessions = [u for u in uniprot_accessions if u not in uniprot2seq.keys()]

print (len(uniprot_accessions), 'uniport accessions to query')

batch_size = 1
batched_uniprot_accessions = [uniprot_accessions[i:i + batch_size] for i in range(0, len(uniprot_accessions), batch_size)]

returned_uniprot = []
for b in tqdm(batched_uniprot_accessions):
    query_string = '+OR+'.join(['accession:{}'.format(s) for s in b])
    res = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={query_string}&format=fasta").text
    returned_uniprot.append(res)


with open('20Nov2023uniprot_chkpoint_for_brenda.txt', 'a') as f:
    for line in returned_uniprot:
        f.write(line)

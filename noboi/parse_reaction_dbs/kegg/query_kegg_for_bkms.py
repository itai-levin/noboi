from bs4 import BeautifulSoup
import re
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
import json
from joblib import Parallel, delayed

bkms_reactions = pd.read_csv('../bkms/1Sep2023_bkms-mapped.txt',sep='\t')
kegg_rids = bkms_reactions['Reaction_ID_KEGG'].drop_duplicates().tolist()[:10]


def get_genes_for_kegg_rid (krid):
    kegg_gene_url = 'https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+rn:'
    url = kegg_gene_url+krid
    res = requests.get(url)
    soup = BeautifulSoup(res.text, 'html.parser')
    return [a.get('href') for a in soup.find_all('a') if a.get('href') and 'entry' in a.get('href')]

def parse_kegg_header (kegg_header):
    if '[' in kegg_header:
        pattern = r'>(\w+:\w+)\s+(\w+)\s+(.*)\s+\['
    else:
        pattern = r'>(\w+:\w+)\s+(\w+)\s+(.*)'
    # Use re.match to apply the pattern to the input string
    match = re.match(pattern, kegg_header)

    if match:
        # Extract matched groups
        gene_id = match.group(1)
        id_value = match.group(2)
        common_name = match.group(3)

        # Create the desired output
        result = {
            'gene_id': gene_id,
            'id': id_value,
            'common_name': common_name,
        }
        ec_pattern = r'\[EC:(.*)\]'
        ec_match = re.findall(ec_pattern, kegg_header)
        if len(ec_match):
            ec_numbers = ec_match[0]
            result['ec_numbers'] = ec_numbers.split(' ')
        return result
    else:
        print("No match found.")
        return None
    
def parse_kegg_fasta(kegg_fasta):
    split = kegg_fasta.split('\n')
    header = split[1]
    results = parse_kegg_header(header)
    results['sequence'] = ''.join(split[2:])
    return results

def retrieve_sequence_for_first_gene(list_of_genes):
    list_of_genes = [x for x in list_of_genes if x]
    if list_of_genes:
        try:
            gene_to_query = list_of_genes[0]
            res = requests.get('https://www.genome.jp/entry/-f+-n+a+'+gene_to_query)
            soup = BeautifulSoup(res.text, 'html.parser')
            return parse_kegg_fasta(soup.find_all('pre')[0].text)
        except IndexError:
            return None
        except TypeError:
            print(list_of_genes)
            return None
        

input_string = ">iho:Igni_0595 K14534 4-hydroxybutyryl-CoA dehydratase / vinylacetyl-CoA-Delta-isomerase [EC:4.2.1.120 5.3.3.3]"
parse_kegg_header(input_string)


unique_kegg_rids =(','.join([str(x) for x in kegg_rids])).split(',')
kegg_genes = [get_genes_for_kegg_rid(x) for x in tqdm(unique_kegg_rids)]
kegg_gene_ids = [[y.split('/')[-1] for y in x] for x in kegg_genes]

kegg_rids_to_genes = dict(zip(unique_kegg_rids, kegg_gene_ids))
with open('kegg_rids_to_genes.json','w') as f:
    json.dump(kegg_rids_to_genes, f)


kegg_sequences = [retrieve_sequence_for_first_gene(x) for x in tqdm(list(kegg_rids_to_genes.values()))]

kegg_rids_to_sequences = dict(zip(unique_kegg_rids, kegg_sequences))
with open('kegg_rids_to_sequences.json','w') as f:
    json.dump(kegg_rids_to_sequences, f)
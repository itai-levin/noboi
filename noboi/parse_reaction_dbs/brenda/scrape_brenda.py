import pandas as pd
from tqdm import tqdm
import requests
from bs4 import BeautifulSoup
import re
import json

class BrendaPage():
   
    def __init__(self, url):
        res = requests.get(url)
        self.soup = BeautifulSoup(res.text, 'html.parser')
        self.get_table_headers()
        self.get_substrate_table_info()

    def get_table_headers(self):
        # return a dictionary of table numbers with a least of the headers for that table
        header_divs = self.soup.find_all('div',id=re.compile('tab[0-9]+_head'))
        self.table_headers_dict = {}
        for tab_div in header_divs:
            tab_num = int(re.findall('(?<=tab)[0-9]+(?=_head)', tab_div.get('id'))[0])
            self.table_headers_dict[tab_num] = [x.strip() for x in tab_div.text.split('\n') if len(x.strip())]

    def _get_substrates_table_id(self):
        #return id for SUBSTRATES table
        titles = ['SUBSTRATE','PRODUCT','REACTION DIAGRAM','ORGANISM','UNIPROT']
        substrate_table_id = [x for x in self.table_headers_dict if self.table_headers_dict[x][:5]==titles]
        return substrate_table_id
    
    def get_substrate_table_info(self):
        substrate_table_id = self._get_substrates_table_id()[0]
        hidden_divs = self.soup.find_all('div', id=re.compile('^tab{}'.format(substrate_table_id)))
        res_list = []
        for div in hidden_divs:
            # Extract values from the first, second, and fifth columns
            res_dict = {}
            for i,c in enumerate(self.table_headers_dict[substrate_table_id]):
                cell = div.find('div', class_='cell', id=lambda x: x and x.endswith('c{}'.format(i)))
                if cell:
                    value = cell.text.strip()
                    res_dict[c] = value
            if res_dict.keys():
                res_list.append(res_dict)
        self.substrate_table_info = res_list

with open('brenda_2023_1.json','r') as f:
    brenda_flatfile = json.load(f)

ec_numbers = list(brenda_flatfile['data'].keys())
urls = ['https://www.brenda-enzymes.org/enzyme.php?ecno={}#SUBSTRATE'.format(ec) for ec in ec_numbers]
brenda_pages_results = {}
for i, (ec, url) in tqdm(enumerate(zip(ec_numbers, urls))):
    try:
        brenda_pages_results[ec] = BrendaPage(url).substrate_table_info
    except Exception as e:
        print (e, 'could not retrieve info for ec number {}'.format(ec))
    if i % 100 == 0:
        with open('scraped_brenda_substrate_results.chkpt', 'w') as f:
            json.dump(brenda_pages_results,f)

with open('scraped_brenda_substrate_results.json', 'w') as f:
    json.dump(brenda_pages_results,f)

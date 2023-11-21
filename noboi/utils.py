from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
import re

def neutralize_atoms(smiles):
    #pulled from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return Chem.MolToSmiles(mol)
    
    
def standardize_smiles(smiles):
    if ' ' in smiles or 'R' in smiles:
        return smiles
    try:
        smiles = neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        return smiles
    except:
        if Chem.MolFromSmiles(smiles)==None:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
            fix_c = AllChem.ReactionFromSmarts('[#6-:1]>>[C;+0:1]')
            if m:
                f = fix_c.RunReactants([m])
                if len (f)>0:
                    f = f[0][0]
                    print ('fixed smiles :', smiles,Chem.MolToSmiles(f))
                    return neutralize_atoms(Chem.MolToSmiles(f))
                else:
                    return smiles
            else:
                return smiles
        else:
            return smiles


def clean_name(string: str) -> str:
    new = re.sub("[\s\.,]", "_", string)
    new = re.sub("[\[\]\(\)']", "", new)
    new = re.sub('&rarr;', '->', new)
    new = re.sub('<[a-zA-z]*>|<\/[a-zA-z]*>|;|&|^[Aa]n |^[0-9]* ','', new)
    new = re.sub('\+', 'plus', new)
    new = re.sub('^-', 'minus', new)
    new = re.sub(',', '-', new)
    return new

def flip_reaction (reaction_smiles):
    reactants, products = reaction_smiles.split('>>')
    return '>>'.join([products, reactants])

def show_rxn_list (list_of_reaction_smiles):
    rd_objs = [AllChem.ReactionFromSmarts(r, useSmiles=True) for r in list_of_reaction_smiles]
    for string, obj in zip(list_of_reaction_smiles, rd_objs):
        plt.figure()
        print(string)
        plt.imshow(Chem.Draw.ReactionToImage(obj))
        plt.show()
        
def standardize_reaction_smiles(reaction_smiles):
    try:
        reactants, products = [standardize_smiles(s) for s in reaction_smiles.split('>>')]
        return '>>'.join([reactants, products])
    except:
        return reaction_smiles

def construct_pathway_from_list (list_of_reactions, metadata = None):
    """
    Inputs:
        list_of_reactions (list of string): list of reaction SMILES
        metadata (list of dict): list of metadata for each reaction. The keys in 
            each dictionary will be used to define the fields of the reaction
            nodes to populate.
    
    Returns:
       pathway tree nested dictionary 
    """
    g = nx.DiGraph()
    
    for i, reaction in enumerate(list_of_reactions):
        reaction = reaction.replace('R','*')
        reactants, intermediate, products = reaction.split('>')
        reactants = [standardize_smiles(s) for s in reactants.split('.')]
        products = [standardize_smiles(s) for s in products.split('.')]
        
        reaction = '>'.join(['.'.join(reactants), intermediate, '.'.join(products)])
        #add reaction nodes
        g.add_node(reaction, color='black', rxn=True, chem=False)
        
        for k in metadata[i].keys():
            g.nodes[reaction][k] = metadata[i][k]

        #connect reactant and product nodes to reaction node
        for reactant in reactants:
            if reactant not in g.nodes:
                g.add_node(reactant, color='red', rxn=False, chem=True, smarts = '')
            g.add_edge(reactant, reaction)      
        for product in products:
            if product not in g.nodes:
                g.add_node(product, color='red', rxn=False, chem=True, smarts = '')
            g.add_edge(reaction, product)
    return g
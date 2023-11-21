from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import networkx as nx
import re
import numpy as np
import pdb


def get_prods_reactants (rxn_str_list, f_or_r):
    """
    Returns a set of string for the reactants, products, and reactions using chemical names

    Parameters:
    rxn_str_list list[str]: reactions with chemical names
    f_or_r (char): 'f' or 'r' whether to look at rxns in the fwd or reverse direction

    Outputs:
    rxtnts (set[string]): of reactant chemical names
    prods (set[string]): product chemical names
    rxns (set[string]): reactions
    """
    prods = set([])
    rxtants = set([])
    rxns = set([])
    simp = re.compile('^[0-9] ')
    for rxn_str in rxn_str_list:
        rxns.add(rxn_str)
        split1 = re.split(' \<?\=\>? ', rxn_str)

        if f_or_r == "f":
            reactants = re.split(" \+ |,", split1[0])
            products = re.split(" \+ |,", split1[1])
        elif f_or_r == "r":
            reactants = re.split(" \+ |,", split1[1])
            products = re.split(" \+ |,", split1[0])
        else:
            print ("forward or backward reactions not specified")
            raise

        for x in reactants:
            s = simp.sub( "", x) #remove articles/ stoichiometric coefficients
            rxtants.add(s)
        for x in products:
            s = simp.sub("", x)
            prods.add(s)
    return rxtants, prods, rxns



def make_coOccurence_tab (rxn_str_list, f_or_r):
    """

    Parameters:
        rxn_str_list list[string]: reactions of form 'Chem_1 [+ Chem_2 + ... ]= Chem 3 ...'
        f_or_r (string): 'f' or 'r' indicates whether to look at rxns in the fwd or reverse direction

    Output:
        numpy array (size num_reactants*num_products) where each entry corresponds to the number of
        reactions in which a  pair of molecules co-occur as a reactant-product pair

    notes: raises an exception if any value other than 'f' or 'r' is passes as f_or_r
    """
    all_reactants, all_products, all_reactions = get_prods_reactants(rxn_str_list, f_or_r)
    all_reactants = list(all_reactants)
    all_products = list(all_products)
    all_reactions = list(all_reactions)
    reactant_to_ind = {}
    product_to_ind = {}
    tot_occurences = np.zeros(len(all_reactants))
    simp = re.compile("^[Aa]n? |^[0-9]+ |^n ")
    for i in range(len(all_reactants)):
        reactant_to_ind[all_reactants[i]] = i
    for i in range(len(all_products)):
        product_to_ind[all_products[i]] = i
    tab = np.zeros((len(all_reactants),len(all_products)))

    for rxn_str in all_reactions:
        split1 = re.split(" \<?\=\>? ", rxn_str)
        if f_or_r == "f":
            reactants = re.split(" \+ |,", split1[0])
            products = re.split(" \+ |,", split1[1])
        elif f_or_r == "r":
            reactants = re.split(" \+ |,", split1[1])
            products = re.split(" \+ |,", split1[0])
        else:
            print ("forward or backward reactions not specified")
            raise
        for reactant in reactants:
            r = simp.sub("", reactant) #remove articles/ stoichiometric coefficients
            r_ind = reactant_to_ind[r]
            tot_occurences[r_ind] += 1
            for prod in products:
                p = simp.sub("", prod) #remove articles/ stoichiometric coefficients
                try:
                    p_ind = product_to_ind[p]
                except Exception as e:
                    print (e)
                tab[r_ind, p_ind] += 1


    return tab, tot_occurences, all_reactants, all_products, reactant_to_ind, product_to_ind

def make_cofactor_dict(rxn_str_list, occ_cutoff, frac_cutoff, f_or_r, ignore=[]):
    """
    Returns a dictionary of 'cofactors,' defined by their appearance in over [occ_cutoff]
    reactions and coappearance with another molecule over [frac_cutoff] of the time
    
    Parameters:
        rxn_str_list (list[string]): reactions represented with common chemical names
        occ_cutoff (int): minimum number of appearances
        frac_cutoff (float): minimum fraction of coappearances
        f_or_r (string): 'f' or 'r' whether to look at rxns in the fwd or reverse direction
    
    output: 
        dict of form :{chemical name: coappearing chemical name}
    """
    
    tab, tot_occurences, reactants, products, r_dic, p_dic = make_coOccurence_tab (rxn_str_list, f_or_r)
    for m in ignore:
        tab[r_dic[m], :] = 0
        tab[:, p_dic[m]] = 0
    r = set()
    p = {}
    rel_tab = tab/tot_occurences[:,None]
    
    #indices for reactants that appear more than 30 times
    for i in np.arange(len(tot_occurences))[tot_occurences > occ_cutoff]: 
        chem = reactants[i]
        if np.max(rel_tab [r_dic[chem],:]) > frac_cutoff:
            r.add(reactants[i])
            l = []
            for j in np.arange(len(rel_tab [r_dic[chem],:])) [np.argsort(rel_tab [r_dic[chem],:])][np.sort(rel_tab [r_dic[chem],:])>frac_cutoff]:
                l.append(products[j])
            p[reactants[i]] = l
    return p
                   



# functions used to process mapped SMILES

def remove_molecule(query_mol, smiles):
    """
    Parameters:
        query_mol (rdkit.Mol): molecule to remove from a molecule SMILES
        smiles (str): molecule SMILES from which to remove molecule

    output:
        SMILES with query_mol removed if it matched one of the molecules
    """
    new_mols_list = []
    for smile in smiles.split('.'):
        mol = Chem.MolFromSmiles(smile)
        if mol:
            if not mol.HasSubstructMatch(query_mol) or not query_mol.HasSubstructMatch(mol):
                new_mols_list.append(mol)
        else:
            return smiles
    return '.'.join([Chem.MolToSmiles(mol) for mol in new_mols_list])

def remove_molecule_from_reactants_and_products (query_mol, reaction_smiles, products_only=False, reactants_only=False):
    """
    Parameters:
        query_mol (rdkit.Mol): molecule to remove from a molecule SMILES
        reaction_smiles (str): reaction SMILES from which to remove moleule
        products_only (bool): if True, remove the molecule only from the products side
        reactants_only (bool):if True, remove the molecule only from the reactants side

    Output
    """
    assert not (products_only and reactants_only)
    if '>' not in reaction_smiles:
        return None
    reactants, intermediate, products = reaction_smiles.split('>')
    if products_only:
        return '>'.join([reactants, intermediate, remove_molecule(query_mol, products)])
    elif reactants_only:
        return '>'.join([remove_molecule(query_mol, reactants), intermediate, products])
    return '>'.join([remove_molecule(query_mol, reactants),intermediate, remove_molecule(query_mol, products)])



def parse_chebi_half_reaction (half_reaction_list, chebi_to_smiles_dict):
    smiles_list = []
    for molecule_chebi in half_reaction_list:
        stoichiometric_coeff = re.findall('^[0-9]+ ', molecule_chebi)
        if len(stoichiometric_coeff):
            stoichiometric_coeff = int(stoichiometric_coeff[0])
        else:
            stoichiometric_coeff = 1

        molecule_chebi = re.sub('^[0-9]+ ', '', molecule_chebi)
        smiles_list += stoichiometric_coeff * [chebi_to_smiles_dict[molecule_chebi]]

    return smiles_list


def remove_cofactors (reactants_list, products_list, cof_dict):
    reactants_to_rm = []
    products_to_rm = []
    for r in reactants_list:
        if r in cof_dict:
            cof_p = cof_dict[r]
            if all([p in products_list for p in cof_p]):
                reactants_to_rm.append(r)
                products_to_rm.extend(cof_p)

    if len(reactants_to_rm) < len(reactants_list) and len(products_to_rm) < len(products_list):
        reactants_list = [r for r in reactants_list if r not in reactants_to_rm]
        products_list = [p for p in products_list if p not in products_to_rm]


    return reactants_list, products_list

def parse_chebi_reaction(reaction, chebi_to_smiles_dict, cof_dict=None, remove_cofs=False):
    reaction_symbol = re.findall(' \<?\=\>? ', reaction)[0]

    chebi_reactants_string, chebi_products_string = reaction.split(reaction_symbol)

    chebi_reactants_list = re.split(' \+ |,', chebi_reactants_string)
    chebi_products_list = re.split(' \+ |,', chebi_products_string)

    if remove_cofs:
        chebi_reactants_list,chebi_products_list=remove_cofactors(chebi_reactants_list,chebi_products_list,cof_dict)
    try:
        reactant_smiles = parse_chebi_half_reaction(chebi_reactants_list, chebi_to_smiles_dict)
    except KeyError as e:
        print ('Could not parse reactants:', e)
        return ''

    try:
        product_smiles = parse_chebi_half_reaction(chebi_products_list, chebi_to_smiles_dict)
    except KeyError as e:
        print ('Could not parse products:', e)
        return ''

    return '.'.join(reactant_smiles)+ '>>' + '.'.join(product_smiles)


def mergeDict(dic1, dic2):
    """
    Merge dictionaries and keep values of common keys in list
    """
    dic = {}
    int_keys = set(dic1).intersection(set(dic2))
    dic3 = {key : list(set(dic1[key]).union(dic2[key])) for key in int_keys}
    dic.update(dic1)
    dic.update(dic2)
    dic.update(dic3)
    return dic

def flip_reaction (reaction_smiles):
    reactants, products = reaction_smiles.split('>>')
    return '>>'.join([products, reactants])

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
                    try:
                        return neutralize_atoms(Chem.MolToSmiles(f))
                    except:
                        return Chem.MolToSmiles(f)
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
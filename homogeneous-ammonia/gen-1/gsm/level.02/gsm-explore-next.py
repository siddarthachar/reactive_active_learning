import sys
import ase
import os
from ase import Atoms
from ase.io import read, write
sys.path.append('/ihome/kjohnson/ska31/AdversarialAttack-DP/reactive_active_learning/')
from se_gsm_wrapper import minimal_wrapper_se_gsm
from deepmd.calculator import DP

import numpy as np
from ase.calculators.orca import ORCA
from ase.calculators.emt import EMT
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from collections import defaultdict
from itertools import combinations

import os
import numpy as np
from ase.calculators.orca import ORCA
from ase.calculators.emt import EMT
from rdkit import Chem
from rdkit.Chem import AllChem
# MyOrcaProfile = OrcaProfile(["/ihome/crc/install/orca/3.0.3/orca"])
# calc = ORCA(profile=MyOrcaProfile)
from ase import Atoms
from ase.io import read, write
import networkx as nx
# import ase.calculators.orca as orc

try:
    from ase import Atoms
    import ase.io
    from ase.io import read, write
except ModuleNotFoundError:
    pass

from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT
from ase.calculators.amber import Amber
from pyGSM.coordinate_systems.delocalized_coordinates import DelocalizedInternalCoordinates
from pyGSM.coordinate_systems.primitive_internals import PrimitiveInternalCoordinates
from pyGSM.coordinate_systems.topology import Topology
from pyGSM.growing_string_methods import DE_GSM
from pyGSM.growing_string_methods import SE_GSM

from pyGSM.level_of_theories.ase import ASELoT
from pyGSM.optimizers.eigenvector_follow import eigenvector_follow
from pyGSM.optimizers.lbfgs import lbfgs
from pyGSM.potential_energy_surfaces import PES
from pyGSM.utilities import nifty
from pyGSM.utilities.elements import ElementData
from pyGSM.molecule import Molecule
from pyGSM.coordinate_systems.slots import Distance, Angle, Dihedral, LinearAngle, MultiAngle
from collections import defaultdict
from itertools import combinations
import pickle
from ase.geometry.analysis import Analysis
# from asdsad import react2prod

reference_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/ral-intermediates.02/gen-1/gsm/level.01/'
prod_dict_path = reference_path + 'prod_0_dicts.pkl'
dp_path='/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/ral-intermediates.02/gen-1/train/dp0/nh.pb'
working_path='/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/ral-intermediates.02/gen-1/gsm/level.02'
total_reactions=6

def adj_matrix2list(adjacency_matrix, remove_repeats=True):
    adjacency_dict = {}
    for i in range(adjacency_matrix.shape[0]):
        neighbors = []
        for j in range(adjacency_matrix.shape[1]):
            if adjacency_matrix[i][j] == 1 and i!=j:
                neighbors.append(j)
        adjacency_dict[i] = neighbors
    return adjacency_dict

def get_mol_symb_list(mol):
    if len([mol]) == 1:
        mol = [mol]
    atom_symbols_list = [atom.GetSymbol() for r in mol for atom in r.GetAtoms()]
    return atom_symbols_list

def get_react_symbols_dict(mol):
    atom_symbols = {}
    atom_symbols_list = [atom.GetSymbol() for r in mol for atom in r.GetAtoms()]
    for i, s in enumerate(atom_symbols_list):
        atom_symbols[i]=s
    return atom_symbols

def get_react_symbols_dict2(mol):
    atom_symb = defaultdict(list)
    atom_symbols_list = [atom.GetSymbol() for r in mol for atom in r.GetAtoms()]
#     print(atom_symbols_list)
    for i, symb in enumerate(atom_symbols_list):
        atom_symb[symb].append(i)
    return atom_symb


def make_reactant_adj_list(reactant_smiles):
    # converts reactant smiles to rdkit mol format
    rec_mol = [Chem.rdmolops.AddHs(Chem.MolFromSmiles(rx)) for rx in reactant_smiles]
    
    # uses rdkit mol format to make adjacency matrix of reach reactant
    adjacency_matrix = [Chem.GetAdjacencyMatrix(molx, useBO=False) for molx in rec_mol] # rdkit does not account for missing Hs while making adj matrix
    
    # makes a list of number of atoms in each reactant
    len_list = [len(x)for x in adjacency_matrix]
    
    # finds the size of the combined adj matrix and makes a big zero matrix
    size_big_react_matrix=sum(len_list)
    big_adj = np.zeros([size_big_react_matrix,size_big_react_matrix])
    
    # combines all reactant adj matrices into the big matrix 
    for i, adj_i in enumerate(adjacency_matrix):
        left_len = sum(len_list[:i]) # len of all atom to the left of current reactant
        right_len = sum(len_list[i+1:])
        pad_matrix = np.pad(adj_i,(left_len,right_len))
        big_adj += pad_matrix
        
    adj_list = adj_matrix2list(big_adj) # adj list of the big reactant matrix
    atom_symbols = get_react_symbols_dict2(rec_mol)
    
    return adj_list, atom_symbols

def make_prod_adj_list(product_smiles):
    prod_mol = [Chem.rdmolops.AddHs(Chem.MolFromSmiles(rx)) for rx in product_smiles]
    adjacency_matrix = [Chem.GetAdjacencyMatrix(molx, useBO=False) for molx in prod_mol]
    adj_list = [adj_matrix2list(x) for x in adjacency_matrix]
    return adj_list

def reindex_dicts(react, prod, react_adj, prod_adj):
    # function that compares the indices of reactant and prod and makes sure 
    # that the prod has indices which are only in the reactant. 
    replacement_dict = {}
    for elem in prod.keys(): # iterat
#         print(f'{elem}: ',react[elem], prod[elem])
        tochange = [x for x in prod[elem] if x not in react[elem]]
        if len(tochange):
            new_prod_list = react[elem][:len(prod[elem])]
            append_replacement_dict = {old:new for old,new in zip(prod[elem],new_prod_list)}
            replacement_dict.update(append_replacement_dict)
            prod[elem]=new_prod_list

    if replacement_dict: # if there is a replacement to be done. else just return as is. 
        old_prod_adj = prod_adj.copy()
        new_dict = {}
        for key, value in prod_adj.items():
            # Replace the key if it exists in the mapping dictionary, or keep it unchanged
            new_key = replacement_dict.get(key, key)
            # Replace values in the list using a list comprehension
            new_value = [replacement_dict.get(item, item) for item in value]
            # Add the updated key-value pair to the new dictionary
            new_dict[new_key] = new_value
        prod_adj = new_dict
    return prod, prod_adj

def generate_add_break_commands(react_adj,prod_adj):
    command_list = []
    for i in prod_adj.keys():
        command_list += [("ADD",*sorted([i+1,j+1])) # +1 to maintain 1 indexing needed for SE-GSM
                         for j in prod_adj[i] 
                         if j not in react_adj[i] 
                         and ("ADD",*sorted([i+1,j+1])) not in command_list]

        command_list += [("BREAK",*sorted([i+1,j+1]))
                         for j in react_adj[i]
                         if j not in prod_adj[i]
                         and ("BREAK",*sorted([i+1,j+1])) not in command_list]
    return(command_list)

# loading the products from level.01
with open(prod_dict_path, 'rb') as f:
    prod_0 = pickle.load(f)

# Uses the products from SE-GSM from the previous level (eg. level 01) and compares it with
# the products from rex_gen. New gsm commands are generated based on this comparison. 
gsm_commands_list=[]
sampled = []
ase_prods_list = []
for i in range(0,total_reactions): # list of all rxns
    gsm_commands=[]
    sampled_rxn = []
    prod=prod_0[i]
    files = []
    ase_prods = []
    rxn_path = reference_path+f'/rxn_{i}'
#     print(os.listdir(rxn_path))
    for j in os.listdir(rxn_path):
        if j.startswith(f'dc_comb'):
            files.append(rxn_path+'/'+j)            
    get_atype = 1
    for k in range(len(files)): # iterate over each dc_comb
        feature = []
        if f'opt_converged_000_ase.xyz' in os.listdir(files[k]):
            # reading only from product for each gsm run
            pos=read(f'{files[k]}/opt_converged_000_ase.xyz',index=-1) 
            sampled_rxn += [files[k]]
            # saving the ase Atoms obj for this product
            ase_prods += [pos]
            # using Analysis to get adjacency matrix of product
            a = Analysis(pos)
            mat = a.adjacency_matrix[0].toarray()
            # converting mat to symmetrix matrix (somehow not before in the prev. step)
            symmetric = mat + mat.T - np.diag(np.diag(mat))
            react=adj_matrix2list(symmetric)
            # comparing the prod from react to get gsm commands
            gsm=generate_add_break_commands(react,prod)
            gsm_commands += [gsm]
    ase_prods_list += [ase_prods]
    sampled += [sampled_rxn]
    gsm_commands_list += [gsm_commands]


def generate_gsm_combinations(data, numadd, numbreak):
#     data = gsm_commands[0]
    add_elements = [item for item in data if item[0] == 'ADD']
    break_elements = [item for item in data if item[0] == 'BREAK']

    # Generate combinations
    combinations_list = []
    for add_comb in combinations(add_elements, numadd):
        for break_comb in combinations(break_elements, numbreak):
            combination = list(add_comb) + list(break_comb)
            combinations_list.append(combination)

    return combinations_list


os.chdir(working_path)
# This initial should be from the previous steps

# the reactants should be sourced from the 
for i, dc_rxn in enumerate(gsm_commands_list):

    os.chdir(working_path)
    if f'rxn_{i}' not in os.listdir():
        os.mkdir(f'rxn_{i}')
    os.chdir(f'rxn_{i}')
    dc_comb_track = -1 # so that it starts from 0 each time
    # need to decompress all dc_combs into a series
    for j, dc in enumerate(dc_rxn):
        if len(dc) > 4:
            # dc_comb_track += 1
            #os.system(f'cp {filepath2} .')
            reactant=ase_prods_list[i][j]
            new_command_3=generate_gsm_combinations(dc,numadd=2,numbreak=1)
            new_command_4=generate_gsm_combinations(dc,numadd=2,numbreak=2)
            new_commands = new_command_3 + new_command_4
            print('****NEW SE-GSM********')
            print(new_commands)
            print(reactant)
            for j in range(len(new_commands)):
                dc_comb_track += 1
    #        for j in range(0,2):
                if f'dc_comb_{dc_comb_track}' not in os.listdir():
                    os.mkdir(f'dc_comb_{dc_comb_track}')
                
                os.chdir(f'dc_comb_{dc_comb_track}')
                print(f"{new_commands[j]}")
                try:
                    minimal_wrapper_se_gsm(reactant, DP(model=dp_path),new_commands[j], num_nodes=20, fixed_reactant=True);  
                    os.chdir('../')
                    
                except Exception as e:
                    print("****got an error*******")
                    print(e)
                    os.chdir('../')
                    continue
                                                     
                
        elif len(dc) <=4 and len(dc) > 0:
            dc_comb_track += 1
            #os.system(f'cp {filepath1} .')
            reactant=ase_prods_list[i][j]
            for k in range(0,1): # just iterating over 1 combination
                if f'dc_comb_{dc_comb_track}' not in os.listdir():
                    os.mkdir(f'dc_comb_{dc_comb_track}')
                os.chdir(f'dc_comb_{dc_comb_track}')
                try:
                    minimal_wrapper_se_gsm(reactant, DP(model=dp_path),dc, num_nodes=20, fixed_reactant=True);  
                    os.chdir('../')
                except Exception as e:
                    print("****got an error*******")
                    print(e)
                    os.chdir('../')
                    continue
        else: # case where there are no commands - the product is achieved
            print('product reached!!')
            continue

            


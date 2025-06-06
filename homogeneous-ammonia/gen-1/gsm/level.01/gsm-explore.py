import sys
import ase
import os
from ase import Atoms
from ase.io import read, write
sys.path.append('/ihome/kjohnson/ska31/AdversarialAttack-DP/reactive_active_learning/')
from se_gsm_wrapper import minimal_wrapper_se_gsm
from deepmd.calculator import DP
from pathlib import Path
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
# from asdsad import react2prod

def adj_matrix2list(adjacency_matrix):
    adjacency_dict = {}
    for i in range(adjacency_matrix.shape[0]):
        neighbors = []
        for j in range(adjacency_matrix.shape[1]):
            if adjacency_matrix[i][j] == 1:
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


reactant_smiles = ['N#N','[H][H]']
add_reactant_smiles = ['N#N','[H][H]'] # if the reactants do not satisfy the products, then we can use these set of reactants to statisfy the stoichiometry.
add_reactant_mol = [Chem.rdmolops.AddHs(Chem.MolFromSmiles(rx)) for rx in add_reactant_smiles]
react_adj_list, react_sym = make_reactant_adj_list(reactant_smiles)
print('Adjacency list of the reactants: \t\t ',react_adj_list)
print('Chemcial symbols and indices of the reactants:   ',react_sym)

product_smiles = ['[NH3]','NN','N=N','[NH+]#[NH+]','N=[NH2+]','[H]']
prod_mol = [Chem.rdmolops.AddHs(Chem.MolFromSmiles(rx)) for rx in product_smiles]

make_prod_adj_list([product_smiles[0]])

# mol = prod_mol[0]
# prod_ind= get_react_symbols_dict2([mol])
gsm_commands = []
for cnt,p in enumerate(prod_mol): # iterate over each product
    prod_dict = get_react_symbols_dict2([p])
    react_sym_i = react_sym.copy()
    react_smiles_i = reactant_smiles.copy()    
    # checking if the elements are conserved:
    elem_checker = [1 if elem in react_sym.keys() else 0 for elem in prod_dict.keys()]
    if 0 not in elem_checker:
        for elem in prod_dict.keys():
            if len(prod_dict[elem]) > len(react_sym[elem]):
                # suggest new reactants to satisfy the defecit and maintain stoichiometry
                # make sure that the add and break command lists are suggested in this section. 
                # append the reactants that can be added to the initial set of reactants and 
                # move the code below to get the command list. 
                
                excess = len(prod_dict[elem])-len(react_sym[elem])
#                 print(excess, elem)
#                 more_react = reactant_smiles.copy()
                
                for add_mol, add_smiles in zip(add_reactant_mol, add_reactant_smiles):
                    # iterate over each additional reactants
                    sym_list = get_mol_symb_list(add_mol) # makes a list for all elems
                    if elem in sym_list:
                        elem_in_add_mol = [x for x in sym_list if x == elem]
                        more_react_elem_list = []
                        while True:
                            # TODO: Try to randomize what reactants to add to the system
                            if len(more_react_elem_list) < excess: 
                                more_react_elem_list += elem_in_add_mol
                                react_smiles_i += [add_smiles]
                            else:
                                break
        print(react_smiles_i, " --> ",product_smiles[cnt])  # note TODO: no more_react for cases with no change - change code to make it consistent                       
        react_adj_list_i, react_elem_ind_i = make_reactant_adj_list(react_smiles_i)
        print('Reactant adjacency: ', react_adj_list_i)
        prod_adj_i = adj_matrix2list(Chem.GetAdjacencyMatrix(p, useBO=False))
        prod_dict, prod_adj_i = reindex_dicts(react_elem_ind_i, prod_dict, react_adj_list_i, prod_adj_i)
        print("Product adjacency: ",prod_adj_i)
        print('Reactant dict: ',react_elem_ind_i)
        commands = generate_add_break_commands(react_adj_list_i, prod_adj_i)
        gsm_commands += [commands]
        print('Commands: ',commands)
        
        # 



# Cases where there are inconsistent number of elements. 
    elif len(prod_dict.keys()) > len(react_sym.keys()): # prod has more elements 
        print('more elements in prod')
        print(prod_dict.keys()-react_sym.keys()) # elements that are in excess
    else:
        print('more elements in react')
        print(react_sym.keys()-prod_dict.keys())
    print('---- \n')
#     break    


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
curr=Path().absolute()
main_path=str(curr.parent.parent)
pick=np.random.randint(4)
dp_path=f'{main_path}/train/dp{pick}/nh.pb'
working_path=os.getcwd()
os.chdir(working_path)
# This initial should be from the previous steps
filepath1 = f"{working_path}/n2-h2.xyz"
reactant1=read(filepath1)
filepath2 = f"{working_path}/n2-h2-h2.xyz"
reactant2=read(filepath2)

for i, dc in enumerate(gsm_commands):
    os.chdir(working_path)
    if f'rxn_{i}' not in os.listdir():
        os.mkdir(f'rxn_{i}')
    os.chdir(f'rxn_{i}')
    if len(dc) > 4:
        os.system(f'cp {filepath2} .')
        reactant=reactant2
        new_command_3=generate_gsm_combinations(dc,numadd=1,numbreak=2)
        new_command_4=generate_gsm_combinations(dc,numadd=2,numbreak=2)
        new_commands = new_command_3 + new_command_4
        print('****NEW SE-GSM********')
        print(new_commands)
        print(reactant)
        for j in range(len(new_commands)):
#        for j in range(0,2):
            if f'dc_comb_{j}' not in os.listdir():
                os.mkdir(f'dc_comb_{j}')
            
            os.chdir(f'dc_comb_{j}')
            print(f"{new_commands[j]}")
            try:
                minimal_wrapper_se_gsm(reactant, DP(model=dp_path),new_commands[j], num_nodes=20, fixed_reactant=True);  
                os.chdir('../')
                
            except Exception as e:
                print("****got an error*******")
                print(e)
                os.chdir('../')
                continue
                                                 
            
    else:
        os.system(f'cp {filepath1} .')
        reactant=reactant1
        for k in range(0,1):
            if f'dc_comb_{k}' not in os.listdir():
                os.mkdir(f'dc_comb_{k}')
            os.chdir(f'dc_comb_{k}')
            try:
                minimal_wrapper_se_gsm(reactant, DP(model=dp_path),dc, num_nodes=20, fixed_reactant=True);  
                os.chdir('../')
            except Exception as e:
                print("****got an error*******")
                print(e)
                os.chdir('../')
                continue

            


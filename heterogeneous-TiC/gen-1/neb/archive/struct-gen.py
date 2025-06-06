from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
from collections import Counter
from ase import Atom
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from deepmd.infer import DeepPot
import os
# from mp_api.client import MPRester
# api_key='VJu07llihiS1dRqq8MPm1OT4p86jWOlc'
from catkit.build import molecule
from catkit.gen.adsorption import Builder
import dpdata
import json
from sklearn.metrics import mean_squared_error
import itertools
import random
from catkit.gen.pathways import ReactionNetwork
from catkit.gen import molecules
from catkit.gen.adsorption import Builder
import networkx as nx
from ase.utils import formula_hill
from ase.db import connect
from ase.optimize import BFGS
from deepmd.calculator import DP
from itertools import combinations_with_replacement
from scipy.spatial.distance import cdist
from ase import Atoms
from scipy.spatial.distance import cdist
from ase.geometry import get_distances

current_gen = 1
# change bulk path as needed. Taking the TiC here. 
bulk_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/bulk/structures/TiC.poscar'
reactant_adsorbate = molecule('C2H6')[0] # change this as needed. This will be stuck on top of surface. Picking the first conformation for now

# NOTE: Remember to change the flag 0. Search for it in the script.

# Adding reactant on surface ###############################

# Create surface with C2H6 adsorbed
bulk_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/bulk/structures/TiC.poscar'
atoms = read(bulk_path)
atoms = atoms.repeat([3,3,1])
gen = SlabGenerator(
    atoms,
    miller_index=(1, 0, 0),
    layers=8,
    fixed=2,
    vacuum=7.5)
slab = gen.get_slab().repeat([3,3,1])
builder = Builder(slab)
# len(builder.get_adsorption_edges()) # use only the length of edges


# generating slab
react_slab = builder.add_adsorbate(molecule(reactant_adsorbate), index = -1, bonds=[0])

# relaxing reactant on slab using DP
dp_path = f'../train/'
for a in react_slab:
    pick = np.random.randint(4)
    a.calc=DP(model=f'{dp_path}/dp{pick}/graph.pb')

for i, a in enumerate(prod_slab):
    dyn = BFGS(a)
    dyn.run(fmax=1e-6)
    write(f'init.struct/ReactantOnSlab-{i}.poscar',a)

# Product enumeration ##################################

## Defining functions

def intermediate_enum(db_name, element_pool):
    '''
    Function that takes in reactant in .db format and enumerates all possible products.
    '''
    with ReactionNetwork(db_name=db_name) as rn:
        # Run a molecule search
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6}, # change this as needed. Maybe it can just be element_pool
            multiple_bond_search=False)

    with ReactionNetwork(db_name=db_name) as rn:
        # Substitution pathway search is expensive!
        rn.path_search(
            reconfiguration=False,
            substitution=False)

    #     rn.plot_reaction_network(file_name='reaction-network.png')
    with ReactionNetwork(db_name=db_name) as rn:
        molecules_list = rn.load_molecules()
        pathways_list = rn.load_pathways()

    all_symbols = []
    for k, v in molecules_list.items():
        all_symbols += [v.get_chemical_formula()]

    unique_symbols = list(set(all_symbols)) # remove repeat symbols
    tops = [molecules.get_topologies(x) for x in unique_symbols] # get all the topologies
    return unique_symbols, tops

def find_combinations(target, counters):
    total_atoms = sum(target.values())
    combos_list = []
    for r in range(1, total_atoms + 1):
        for combo in combinations_with_replacement(counters, r):
#             print(combo)
            current = sum(combo, Counter())
            if current == target:
                print(list(combo)) # convert to string here and return list of list of strings.
                combos_list += [list(combo)]

    return combos_list

def get_formula_from_counter(x, elem_order):
    p=''
    for e in elem_order:
        cnt = x[e]
        if cnt == 1:
            toadd = str(e)
        elif cnt == 0:
            continue 
        else:
            toadd = f'{e}{cnt}'
        p += toadd
    return p

# We need to sample all possible product degredation reactions with all the intermediates we have
# need to start with C2H6 and CH4
db_name = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/molecule/ethane.db' # this was generated using ASE. Just save as .db.
sym, tops = intermediate_enum(db_name, {'C':2, 'H':6})
# all the intermediates for C2H6 - we need the symbol

reactant_adsorbate_s = reactant_adsorbate.get_chemical_symbols() # chemical symbols
adsorbate_sc = Counter(reactant_adsorbate_s) # 
int_s = [Counter(t[0].get_chemical_symbols()) for t in tops] # get all the dict for each symbol (interm.)

combos_list = find_combinations(target_counter, int_s) # generates combinations of products to match reactants

# generating formula from combinations
mols = []
for x in combos_list:
    mols.append([get_formula_from_counter(y, ['C','H']) for y in x])
print(mols)

# Reloading empty surface
# Adding indiv molecules
atoms = read(bulk_path)
atoms = atoms.repeat([3,3,1])
gen = SlabGenerator(
    atoms,
    miller_index=(1, 0, 0),
    layers=8,
    fixed=2,
    vacuum=7.5)
slab = gen.get_slab().repeat([3,3,1])
builder = Builder(slab)
# len(builder.get_adsorption_edges()) # use only the length of edges

def set_tags_random(x):
    tags = x.get_tags()
    pick = random.randrange(len(tags))
    tags[pick] = -1
    x.set_tags(tags)
    return x

save_adsorbates = 'init.struct'
if 'add_adsorbates' not in os.listdir(save_adsorbates):
    os.mkdir(f'{save_adsorbates}/add_adsorbates')

for i, mol in enumerate(mols):
    mol_obj = [molecule(x)[0] for x in mol]
    tags = [x.get_tags() for x in mol_obj]
    # randomly set one of the atoms to be -1 and see if it can be adsorbed on the site
    mols_update = [set_tags_random(x) for x in mol_obj]
    builder = Builder(slab) # make new slab
    ads_slab = builder.add_adsorbates(mols_update, indices = list(np.arange(len(mols_update))))
    write(f'{save_adsorbates}/add_adsorbates/ProdOnSurface-{str(i).zfill(3)}.poscar',ads_slab, sort = False)

# Reindexing products based on reactants  ##################

if 'add_adsorbates_mol_only_reorder' not in os.listdir(save_adsorbates):
    os.mkdir(f'{save_adsorbates}/add_adsorbates_mol_only_reorder')

final_add_ads = f'{save_adsorbates}/add_adsorbates/'
react = f'{save_adsorbates}/ReactantOnSlab-0.poscar' #flag 0
prods_i = [x for x in os.listdir(f'{save_adsorbates}/add_adsorbates')] 
reoreder_save = f'{save_adsorbates}/add_adsorbates_mol_only_reorder'
prods = [f'{save_adsorbates}/add_adsorbates/{x}' for x in os.listdir(f'{save_adsorbates}/add_adsorbates')]
r = read(react) # loading reactant 
# convert to atoms object (even for products) and use the rearranging code.
# we can do an element wise screening 

# Functions
def create_atom_from_pos_symb(mol_symb, mol_pos):
#     mol_symb = r.get_chemical_symbols()[-8:]
#     mol_pos = r.get_positions()[-8:]
    symbs = list(dict.fromkeys(mol_symb))
    formula = ''.join([f'{symbs[i]}{len([x for x in mol_symb if x == symbs[i]])}' for i in range(len(symbs))])
    mol = Atoms(formula)
    mol.set_chemical_symbols(mol_symb)
    mol.set_positions(mol_pos)
    return mol

def reorder_atoms(a_atoms, b_atoms):
    # Extract atomic symbols and positions
    a_symbols = a_atoms.get_chemical_symbols()
    b_symbols = b_atoms.get_chemical_symbols()
    a_positions = a_atoms.get_positions()
    b_positions = b_atoms.get_positions()

    # Calculate pairwise distances between atoms in A and B
    distances = cdist(a_positions, b_positions)

    # Create a mapping from indices of atoms in B to indices of atoms in A
    atom_mapping = []
    used_b_indices = set()
    
    for i in range(len(a_atoms)):
        # Find the closest unmatched atom in B to the ith atom in A
        min_dist_idx_a, min_dist_idx_b = np.unravel_index(np.argmin(distances), distances.shape)
        
        # Ensure that the atoms have the same species
        while b_symbols[min_dist_idx_b] != a_symbols[min_dist_idx_a] or min_dist_idx_b in used_b_indices:
            distances[min_dist_idx_a, min_dist_idx_b] = np.inf
            min_dist_idx_a, min_dist_idx_b = np.unravel_index(np.argmin(distances), distances.shape)
        
        atom_mapping.append((min_dist_idx_a, min_dist_idx_b))
        used_b_indices.add(min_dist_idx_b)
        
        distances[min_dist_idx_a, :] = np.inf
        distances[:, min_dist_idx_b] = np.inf

    # Sort atom_mapping based on the index in A
    atom_mapping.sort()

    # Rearrange atoms in B based on the mapping
    b_reordered_positions = [b_positions[idx_b] for _, idx_b in atom_mapping]
    b_reordered_symbols = [b_symbols[idx_b] for _, idx_b in atom_mapping]
    
    b_reordered = Atoms(symbols=b_reordered_symbols, positions=b_reordered_positions)
    
    return b_reordered


## final code to reorder all prods. 
# Used the final 8 atoms here because of reactant (ethane).
for i, ppath in enumerate(prods):
    p = read(ppath)
    org_chem_sym = p.get_chemical_symbols()
    org_pos = p.get_positions()
    p_mol = create_atom_from_pos_symb(org_chem_sym[-8:], org_pos[-8:])
    p_ordered = reorder_atoms(react_mol, p_mol)
    org_chem_sym[-8:] = p_ordered.get_chemical_symbols()
    org_pos[-8:] = p_ordered.get_positions()
    p.set_positions(org_pos)
    p.set_chemical_symbols(org_chem_sym)
    write(f'{reoreder_save}/ProdOnSurface-{str(i).zfill(3)}.poscar', p)

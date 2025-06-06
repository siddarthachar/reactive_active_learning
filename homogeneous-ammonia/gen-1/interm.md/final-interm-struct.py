import numpy as np
import os
from ase.io import read, write
import mddatasetbuilder
from mddatasetbuilder.detect import DetectDump
kB = 8.6173303e-5 #eV/K
from deepmd.calculator import DP
from ase.optimize import BFGS, FIRE
from rdkit import Chem
import MDAnalysis as mda
import mdapackmol
import dpdata

# things to take care of in this code: 
# 1. box size and number of molecules to pack. 
# 2. Make sure to have your molecules stored as indiv xyz files in ../../init.mol
# 3. Change the "desired order" in the second last time when fixing the ordering in POSCAR.

def bonds_nopbc(atoms):
    bonds = DetectDump._crd2bond(atoms, False)
    levels = DetectDump._crd2bond(atoms, True)
    return bonds, levels

def get_multiplicity(valence_e, atm):
    # valence_e = {'N':5, 'H':1}
    sym = atm.get_chemical_symbols()
    elem_track = {i:valence_e[s] for i, s in enumerate(sym)}
    bonds, btype = bonds_nopbc(atm)
    print(btype)
    for i, bs in enumerate(zip(bonds, btype)):
        b, bt = bs
        paired_e = sum(bt)
        if sym[i] == 'O' and len(bt) == 1: # treating O2 as triplet state (for reactions)
            if bt[0] == 2 and sym[b[0]] == 'O': # making sure that the O is double bonded to another O
                paired_e = 1
        elem_track[i] -= paired_e # summing the total bonds
    total_unpaired = sum(elem_track.values())
    multiplicity = total_unpaired + 1
    return multiplicity, total_unpaired

def relax_atoms(atoms, dp_path, pb_name='nh.pb', logpath='rlx.out'):
    last.calc = DP(f'{dp_path}/dp{np.random.randint(4)}/{pb_name}')
    relax = BFGS(last, logfile = logpath)
    relax.run(fmax=0.01)
    return relax.atoms

def poscar2lmp(inp, out):
    data = dpdata.System(inp, fmt = 'vasp/poscar')
    data.to('lammps/lmp', out, frame_idx=0)
    
def maintain_order_poscar(atm, inp, desired_order):
    # inp needs the path.

    # atm = read(inp)
    # np.unique(atm.get_chemical_symbols())
    tot_atms_list = [len([x for x in atm.get_chemical_symbols() if x == do]) for do in desired_order] # [96, 48]
    desired_chemf = ''.join([str(x) + str(y) for x, y in zip(desired_order, tot_atms_list)]) # desired chem formula
    og_pos = atm.get_positions()
    if desired_chemf != atm.get_chemical_formula():
        print('pain --------------')
    #     print('pain')
        atm_2 = atm.copy()
        new_symb = list(itertools.chain.from_iterable(itertools.repeat(x, y) for x,y in zip(desired_order, tot_atms_list)))
        atm_2.set_chemical_symbols(new_symb)
        new_pos = np.vstack([np.array([x for x, y in zip(atm.get_positions(),atm.get_chemical_symbols())  if y == do]) for do in desired_order])
        atm_2.set_positions(new_pos)
    #     print(atm_2)
        write(f'{inp}',atm_2)
        print(f'wrote to {inp}')


# script to iterate over all gsm reactions and find the barrier and
#the product (with the product's SMILES and )
# interm_dict = {}
path = '../gsm/'
dp_path = '../train'
# interm = []
rlx_interm = []
T = 300 # K
kBT = 2*kB*T # using Boltzmann factor to account for the likelihood or stability of the product. 
init_mols_path = '../../init.mol' # for packmoling

if 'struct.interm' not in os.listdir():
    os.mkdir('struct.interm')
if 'struct.packed' not in os.listdir():
    os.mkdir('struct.packed')

struct_count = 0
for level in sorted(os.listdir(path)):
    li = f'{path}/{level}/'
    for rxn in sorted([x for x in os.listdir(li) if 'rxn' in x]):
        print(rxn)
        rxn_path = f'{li}/{rxn}'
        for dc in sorted([x for x in os.listdir(rxn_path) if 'dc_comb' in x]):
            dc_path = f'{rxn_path}/{dc}/'
            if 'opt_converged_000_ase.xyz' in os.listdir(dc_path):
                ob = read(f'{dc_path}/opt_converged_000_ase.xyz', index=':')
                react = ob[0]
                react_bonds = bonds_nopbc(react) # bond adj list of reactant
                max_en = max([x.get_potential_energy() for x in ob])
                if max_en-ob[-1].get_potential_energy() > kBT:

                    last = ob[-1].copy()
                    initial_bonds = bonds_nopbc(last)
#                     print(react_bonds[0] == initial_bonds[0])
                    if react_bonds[0] != initial_bonds[0]: # if the initial and final are NOT the same. 
                        
                        # interm += [last] # saveing the interm (pre-relaxed)                        
                        relaxed = relax_atoms(last, dp_path, 'nh.pb') # relaxing the final image from SE-GSM
                        # change nh.pb to the protobuf format that you have. 
                        relaxed_bonds = bonds_nopbc(relaxed) # getting the bond adj list
                        # deleting atoms that are radicals and saving them. 
                        del relaxed[[i for i, x in enumerate(relaxed_bonds[0]) if not len(x)]] # atoms to be deleted - if they are radicals post relax
                        rlx_interm += [relaxed]
                        _,unp = get_multiplicity({'N':5,'H':1},relaxed)
                        isnotradical = unp%2 == 0 # checking if it's a radical based on even number of unpaired electrons or not. 
                        if isnotradical :
                            
                            write(f'struct.interm/{str(struct_count).zfill(3)}.xyz', relaxed)
                            # put additional conditions for valence check # added to the top. 

                            
                            # initial set of molecules for the box
                            init_mols = [mda.Universe(f'{init_mols_path}/{x}') for x in os.listdir(init_mols_path) if 'xyz' in x]
                            init_mols += [mda.Universe(f'struct.interm/{str(struct_count).zfill(3)}.xyz')]
                            system = mdapackmol.packmol(
                                [mdapackmol.PackmolStructure(
                                    x, number=10,
                                    instructions=['inside box 0. 0. 0. 20. 20. 20.']) for x in init_mols]
                            )

                            system.atoms.write(f'struct.packed/{str(struct_count).zfill(3)}.xyz')
                            packed_ase = read(f'struct.packed/{str(struct_count).zfill(3)}.xyz') 
                            packed_ase.set_cell([20,20,20])
                            write(f'struct.packed/{str(struct_count).zfill(3)}.poscar', packed_ase, sort=True) # save to poscar
                            maintain_order_poscar(packed_ase, f'struct.packed/{str(struct_count).zfill(3)}.poscar', ['H','N'])
                            poscar2lmp(f'struct.packed/{str(struct_count).zfill(3)}.poscar', f'struct.packed/{str(struct_count).zfill(3)}.lmp')
                            struct_count += 1

print('Completed! ;)')

# # Performing the packmoling bit

# # load individual molecule files
# init_mols_path = '../../init.mol'
# init_mols = [mda.Universe(f'{init_mols_path}/{x}') for x in os.listdir(init_mols_path) if 'xyz' in x]
# init_mols += [mda.Universe(f'struct.interm/{str(0).zfill(3)}.xyz')]
# # water = mda.Universe('water.pdb')
# # urea = mda.Universe('urea.pdb')

# # call Packmol with MDAnalysis objects as arguments
# # the 'instructions' allow for any valid Packmol commands
# system = mdapackmol.packmol(
#     [mdapackmol.PackmolStructure(
#         x, number=10,
#         instructions=['inside box 0. 0. 0. 20. 20. 20.']) for x in init_mols]

#     # [mdapackmol.PackmolStructure(
#     #     water, number=1000,
#     #     instructions=['inside box 0. 0. 0. 40. 40. 40.']),
#     #  mdapackmol.PackmolStructure(
#     #     urea, number=400,
#     #     instructions=['inside box 0. 0. 0. 40. 40. 40.'])]
# )

# system.atoms.write('output.xyz')


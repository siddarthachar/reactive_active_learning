from ase.io import read, write 
import numpy as np
import os 
from pathlib import Path

# structpath = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/reaction/test/'
# prods = os.listdir(f'{structpath}/add_adsorbates_mol_only_reorder/')
# savepath = f'{structpath}/neb/final_neb/'

curr=Path().absolute()
main_path=str(curr.parent)
structpath = f'{main_path}/init.struct/'

reactants=[x for x in os.listdir(structpath) if 'ReactantOnSlab' in x]
prods = os.listdir(f'{structpath}/add_adsorbates_mol_only_reorder')
for i, r in enumerate(reactants):
	nebpath=f'{str(i).zfill(2)}'
	os.chdir(nebpath)
	if 'final_neb' not in os.listdir():
		os.mkdir('final_neb')
	for j in range(len(prods)):
		a=read(f'neb.{str(j).zfill(3)}/neb.traj',index='-7:')
		write(f'final_neb/neb_{str(j).zfill(3)}.traj',a)
	os.chdir('../')

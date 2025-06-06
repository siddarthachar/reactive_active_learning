from ase.io import read, write 
import numpy as np
import os 

structpath = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/reaction/test/'
prods = os.listdir(f'{structpath}/add_adsorbates_mol_only_reorder/')
savepath = f'{structpath}/neb/final_neb/'
if 'final_neb' not in os.listdir():
	os.mkdir('final_neb')
for i in range(len(prods)):
	a=read(f'{i}/neb.traj',index='-7:')
	write(f'{savepath}/neb_{str(i).zfill(3)}.traj',a)

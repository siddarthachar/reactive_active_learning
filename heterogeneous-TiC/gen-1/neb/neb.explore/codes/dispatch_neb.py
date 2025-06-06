import os
import numpy as np
from pathlib import Path


curr=Path().absolute()
main_path=str(curr.parent.parent)
structpath = f'{main_path}/init.struct/'
prods = os.listdir(f'{structpath}/add_adsorbates_mol_only_reorder/')
for i in range(len(prods)):
	nebpath=f'neb.{str(i).zfill(3)}'
	if f'{nebpath}' not in os.listdir():
		os.mkdir(f'{nebpath}')
	os.system(f'cp run_neb.py {nebpath}/')
	os.system(f'cp job.slurm {nebpath}/')
	os.chdir(f'{nebpath}')
	os.system('rm neb.traj')
	os.system('rm slurm*')
	os.system('sbatch job.slurm')
	os.chdir('../')

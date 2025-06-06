import os
import numpy as np

structpath = os.getcwd()
prods = os.listdir(f'{structpath}/add_adsorbates_mol_only_reorder/')
for i in range(len(prods)):
	if f'neb.{str(i).zfill(3)}' not in os.listdir():
		os.mkdir(f'neb.{str(i).zfill(3)')
	os.system(f'cp run_neb.py {i}/')
	os.system(f'cp job.slurm {i}/')
	os.chdir(f'{i}')
	os.system('rm neb.traj')
	os.system('rm slurm*')
	os.system('sbatch job.slurm')
	os.chdir('../')

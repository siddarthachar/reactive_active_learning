import os
from pathlib import Path

curr=Path().absolute()
main_path=str(curr.parent)
structpath = f'{main_path}/init.struct/'

reactants=[x for x in os.listdir(structpath) if 'ReactantOnSlab' in x]

for i, r in enumerate(reactants):
	nebpath=f'{str(i).zfill(2)}'
	if nebpath not in os.listdir():
		os.mkdir(nebpath)
	os.system(f'cp codes/dispatch_neb.py {nebpath}/')
	os.system(f'cp codes/run_neb.py {nebpath}/')
	os.system(f'cp codes/job.slurm {nebpath}/')
	os.chdir(nebpath)
	os.system('python dispatch_neb.py')
	os.chdir('../')
import os

md_paths = sorted([x for x in os.listdir() if 'md' in x])

for md in md_paths:
	os.chdir(md)
	os.system('sbatch job.slurm')
	os.chdir('../')

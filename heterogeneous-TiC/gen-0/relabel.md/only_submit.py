import os

wrkdir = os.getcwd()
md_list = [x for x in os.listdir() if 'md' in x]
for md in md_list:
	image_list = [x for x in os.listdir(md) if 'image' in x]
	for im in image_list:
		os.chdir(f'{wrkdir}/{md}/{im}')
		os.system('sbatch job.slurm')
		os.chdir(wrkdir)
	

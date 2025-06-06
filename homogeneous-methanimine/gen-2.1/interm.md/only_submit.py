import os

md_paths = sorted([x for x in os.listdir() if 'md' in x])
sub_path=open('submitted_md.path','w+')
i = 0
submitted = []
for md in md_paths:
	os.chdir(md)
	submitted += [md]
	os.system('sbatch job.slurm')
	os.chdir('../')
	sub_path.write(f'{md}\n')

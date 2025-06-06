import os
import time
sleep_time = 0
wrkdir = os.getcwd()
md_list = [x for x in os.listdir() if 'md' in x and os.path.isdir(x)]
for md in md_list:
	image_list = [x for x in os.listdir(md) if 'image' in x]
	for im in image_list:
		os.chdir(f'{wrkdir}/{md}/{im}')
		os.system('cp /bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/md-al/sample-vasp/job.slurm .')
		if 'report.out' not in os.listdir():
			os.system('sbatch job.slurm')
			print(f'{md}/{im}')
			sleep_time += 1
		elif len(os.popen('grep JobId: report.out').read())==0:
			os.system('sbatch job.slurm')
			sleep_time += 1
			print(f'{md}/{im}')
		os.chdir(wrkdir)
		if sleep_time == 750:
			time.sleep(60*2)
			sleep_time = 0

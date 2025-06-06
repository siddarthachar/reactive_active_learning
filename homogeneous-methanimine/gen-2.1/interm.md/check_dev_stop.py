import os
import numpy as np
import time
import subprocess

# reading job IDs
job_path = os.popen('cat jobs.dat').read().split('\n')[:-1]
job_ids = [int(x.split()[3]) for x in job_path]
md_paths = os.popen('cat submitted_md.path').read().split('\n')[:-1]

all_finished = []

while len(all_finished)<len(md_paths): # stop when all jobs are finished

	for md, jid in zip(md_paths, job_ids):
		# a = np.loadtxt(f'{md}/model_devi.out')[:,4][-1]
		if md not in all_finished:
			try:
				maxforce = np.loadtxt(f'{md}/model_devi.out')[:,4][-1]
				if maxforce > 10:
					os.system(f'scancel -M gpu {jid}')
					print(f'{jid} {md}: Cancelled')
					all_finished += [md]
				else: # cases when maxforce < 10 till end of simulation
					finished = len(os.popen(f'grep wall {md}/log.lammps')) > 0

					if finished:
						all_finished += [md]
					else:
						print(f'{jid} {md}: Still Running')
			except:
				print(f'{md} {jid} Not started')
	time.sleep(10)
	print('Finished: ', all_finished)


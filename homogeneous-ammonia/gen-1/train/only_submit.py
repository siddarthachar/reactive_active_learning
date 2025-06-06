import os
import json
from random import randint, randrange

main_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH'
current_gen=1 # if current gen is 0 then copy input.json from main path
for i in range(4):
	os.chdir(f'dp{i}')
	os.system('sbatch job.slurm') # remember to change the init model in the job slrum
	os.chdir('../')


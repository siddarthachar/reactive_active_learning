import os
import json
from random import randint, randrange

bulk_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/bulk'
surf_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/surface'
molecule_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/molecule' 

for i in range(4):
        if f'dp{i}' not in os.listdir():
                os.mkdir(f'dp{i}')
        os.system(f'cp input.json dp{i}/')
        os.system(f'cp job.slurm dp{i}/job.slurm') # No init model in job.slurm

        path=f'dp{i}/input.json'
        deepmd_data_curr_b = os.popen(f'cat {bulk_path}/init.data/deepmd_data_path.dat').read().split('\n')[:-1]
        deepmd_data_curr_s = os.popen(f'cat {surf_path}/init.data/deepmd_data_path.dat').read().split('\n')[:-1]
        deepmd_data_curr_m = os.popen(f'cat {molecule_path}/init.data/deepmd_data_path.dat').read().split('\n')[:-1]
        with open(path, 'r+') as f:
            data = json.load(f)
            data["model"]["fitting_net"]["seed"] = randint(1e6, 1e7)
            data["model"]["descriptor"]["seed"] = randint(1e6, 1e7)
            data["training"]["seed"] = randint(1e6, 1e7)
            data["training"]["training_data"]["systems"] = deepmd_data_curr_b + deepmd_data_curr_s + deepmd_data_curr_m
            data["training"]["stop_batch"] = 400000
            data["training"]["training_data"]["batch_size"] = "auto"

        os.remove(path)
        with open(path, 'w') as f:
            json.dump(data, f, indent = 4)
        os.chdir(f'dp{i}')
        os.system('sbatch job.slurm') # remember to change the init model in the job slrum
        os.chdir('../')

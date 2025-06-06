import os
import json
from random import randint, randrange

def replace_nth_line(file_path, n, replacement_string):
    # Read the contents of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Check if the specified line number is valid
    if 1 <= n <= len(lines):
        # Replace the nth line with the new string
        lines[n] = replacement_string + '\n'
        
        os.system(f'mv {file_path} {file_path}.bk')
        
        # Write the modified content back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)
        print(f"Successfully replaced line {n} with: {replacement_string}")
    else:
        print(f"Invalid line number: {n}. The file has {len(lines)} lines.")

# main_path = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH'
main_path = '/bgfs/kjohnson/pbs13/DeePMD/reactive_active_learning/methyle_03/new_DP0/ral.self_constistent/'
current_gen=2.1 # if current gen is 0 then copy input.json from main path
for i in range(4):
	if f'dp{i}' not in os.listdir():
		os.mkdir(f'dp{i}')
	if current_gen:
		os.system(f'cp {main_path}/gen-1.10_prime/train/dp{i}/input.json dp{i}/')
	#os.system(f'cp {main_path}/dp.slurm dp{i}/job.slurm')

	path=f'dp{i}/input.json'
	levels = [x for x in sorted(os.listdir(f'../../gen-{int(round(current_gen-0.1,1))}/relabel.gsm/')) if x.split('.')[0] == 'level']
	deepmd_data_curr = []

	for level in levels:
		deepmd_data_curr += os.popen(f'cat ../../gen-{int(round(current_gen-0.1,1))}/relabel.gsm/{level}/deepmd_data_path.dat').read().split('\n')[:-1]

	# including relabeled md configurations
	deepmd_data_curr += os.popen(f'cat ../../gen-{int(round(current_gen-0.1,1))}/relabel.md/deepmd_data_path.dat').read().split('\n')[:-1]

	with open(path, 'r+') as f:
		data = json.load(f)
		app_sys =  data["training"]["training_data"]["systems"] + deepmd_data_curr
		data["model"]["fitting_net"]["seed"] = randint(1e6, 1e7)
		data["model"]["descriptor"]["seed"] = randint(1e6, 1e7)
#		data["model"]["descriptor"]["sel"] = [92, 46]
#		data["model"]["descriptor"]["rcut"] = 6.0
#		data["model"]["descriptor"]["rcut_smth"] = 2.0
#		data["model"]["descriptor"]["neuron"]=[120,120,120]
#		data["model"]["fitting_net"]["neuron"]=[240,240,240]
		data["training"]["seed"] = randint(1e6, 1e7)
		data["training"]["training_data"]["systems"] = app_sys
		data["training"]["stop_batch"] = 1000000
		data["training"]["training_data"]["batch_size"] = "auto"

	os.remove(path)
	with open(path, 'w') as f:
	    json.dump(data, f, indent = 4)
	os.chdir(f'dp{i}')
	os.system(f'cp {main_path}/gen-{int(round(current_gen-0.1,1))}/train/dp{i}/job.slurm .')
	job_init_model_train = f'dp train input.json --init-model {main_path}/gen-{int(round(current_gen-0.1,1))}/train/dp{i}/model.ckpt'
	replace_nth_line('job.slurm', 13, job_init_model_train) # 13th line is where we have the dp train command
#	os.system('sbatch job.slurm') # remember to change the init model in the job slurm
	os.chdir('../')

	




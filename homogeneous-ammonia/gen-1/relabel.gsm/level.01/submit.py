import os

# submit jobs

#main_path='/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/gen-0/relabel/'
main_path = os.getcwd()
data_list=sorted([x for x in os.listdir(f'{main_path}') if x.split('.')[0]=='rxn'])

print(data_list)
for i, d in enumerate(data_list):
    data_path = f'{main_path}/{d}'
    image_list = sorted(os.listdir(data_path))
    #image_list = sorted([x for x in os.listdir(f'{data_path}') if x.split('.')[0]=='image'])
    for j, im in enumerate(image_list):
        image_path = f'{data_path}/{im}'
        os.chdir(image_path)
        if "OUTCAR" not in os.listdir():
            print(image_path)
            os.system('sbatch job.slurm')
        os.chdir(main_path)

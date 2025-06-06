import os
import numpy as np
import ase
from ase.io import read, write
# Need to add two things: 
# a.  Do not include structures where the SCF has not converged in 200 steps. 
# b.  Modify the filtering code to account for this criterion. 
# c.  When total magnetization is non integer. 

# Atomic energy of H (gamma only): -1.11036 eV, mag = 1.0
# Atomic energy of N (gamma only): -3.10869 eV, mag = 3.0

atom_e = {'N':-3.10869, 'H':-1.11036}

deepmd_paths = []
main_path=os.getcwd()
data_list=sorted([x for x in os.listdir(f'{main_path}') if x.split('.')[0]=='rxn'])
# data_list=sorted([x for x in os.listdir(f'{main_path}') if x.split('.')[0]=='data'])
for i, d in enumerate(data_list):
    data_path = f'{main_path}/{d}'
    if 'deepmd_data' in os.listdir(data_path):
        os.system(f'mv {data_path}/deepmd_data {data_path}/deepmd_data.{np.random.randint(1000)}.bk')
    image_list = sorted([x for x in os.listdir(data_path) if x != 'deepmd_data'])
    # image_list = sorted([x for x in os.listdir(data_path) if x.split('.')[0] == 'image'])
    outcars = []
    for j, im in enumerate(image_list):
        image_path = f'{data_path}/{im}'
        config = read(f'{image_path}/POSCAR')
        elem = config.get_chemical_symbols()
        total_atom_sum = sum([atom_e[i] for i in elem])
        # tests
        contcar_test = os.path.getsize(f'{image_path}/CONTCAR') > 0
        scf_test = int(os.popen(f'grep RMM: {image_path}/OSZICAR').read().split('\n')[-2].split()[1]) < 200 # less than 200 scf steps
        mag_read = os.popen(f'grep F= {image_path}/OSZICAR').read().split()[-1]
        mag_int_check = (float(mag_read) - int(float(mag_read))) < 1e-3 # checking if total magnetic moment is within a given threshold or not.
        # mag_int_check = float(os.popen(f'grep F= {image_path}/OSZICAR').read().split()[-1]).is_integer() # checking if total magnetic moment is integer or not.
        atomization_energy_check = float(os.popen(f'grep F= {image_path}/OSZICAR').read().split()[2]) < total_atom_sum # checking is total energy is less than the total sum of atomic energies.
        # print(contcar_test + scf_test + mag_int_check + atomization_energy_check)
        # print(os.path.getsize(f'{image_path}/CONTCAR'))
        # print(int(os.popen(f'grep RMM: {image_path}/OSZICAR').read().split('\n')[-2].split()[1]))
        # print(float(os.popen(f'grep F= {image_path}/OSZICAR').read().split()[-1]))
        # print(float(os.popen(f'grep F= {image_path}/OSZICAR').read().split()[-1]))
        if (contcar_test + scf_test + mag_int_check + atomization_energy_check) == 4: # if all of them are true
            print(image_path)
            a=read(f'{image_path}/OUTCAR')
            outcars += [a]
    box = []
    coord = []
    energy = []
    force = []
    for i, o in enumerate(outcars):
        if i == 0:
            a=o.get_chemical_symbols()
            indexes = np.unique(a, return_index=True)[1]
            type_map = [a[index] for index in sorted(indexes)]
            type_file = [type_map.index(x) for x in a]
        e = o.get_potential_energy()
        if e < 0: # eV # this is redundant now, but whatever. 
            box += [np.array(o.cell).ravel()]
            coord += [o.get_positions().ravel()]
            energy += [e]
            force += [o.get_forces().ravel()]
    energy = np.array(energy)
    force = np.array(force)
    coord = np.array(coord)
    box = np.array(box)
    deepmd_data = f'{data_path}/deepmd_data'
    deepmd_paths += [deepmd_data]
    set_save = deepmd_data + '/set.000'
    if 'deepmd_data' not in os.listdir(data_path):
        os.mkdir(deepmd_data)
    if 'set.000' not in os.listdir(deepmd_data):
        os.mkdir(set_save)
    np.save(set_save+'/energy.npy', energy)
    np.save(set_save+'/force.npy', force)
    np.save(set_save+'/coord.npy', coord)
    np.save(set_save+'/box.npy', box)
    type_map_file = open(deepmd_data+'/type_map.raw','w+')
    for i, elem in enumerate(type_map):
        type_map_file.write(str(elem)+'\n')
    type_map_file.close()
            
    type_file_o = open(deepmd_data+'/type.raw','w+')
    for i, elem in enumerate(type_file):
        type_file_o.write(str(elem)+'\n')
    type_file_o.close()

dp_files = open(f'{main_path}/deepmd_data_path.dat', 'w+')
for i, pp in enumerate(deepmd_paths):
    dp_files.write(str(pp)+'\n')
dp_files.close()

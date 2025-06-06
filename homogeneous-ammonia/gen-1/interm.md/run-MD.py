import os
from ase.io import read, write
import numpy as np


md_path = os.getcwd()
dp_path = '../train/'
for i in range(4):
    os.system(f'ln -s {dp_path}/dp{i}/nh.pb {md_path}/graph.{str(i).zfill(3)}.pb')

# collecting lmp_files
lmp_files = sorted([f'{md_path}/struct.packed/'+x for x in os.listdir('struct.packed') if 'lmp' in x])


ref_lammps = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/ANH/md-al/sample_lammps'
lmp_12 = f'{ref_lammps}/in-12.lammps' # HN

for i, lm in enumerate(lmp_files):
    if f'md.{str(i).zfill(2)}' not in os.listdir(md_path):
        os.mkdir(f'{md_path}/md.{str(i).zfill(2)}')
    md_i = f'{md_path}/md.{str(i).zfill(2)}'
    vsp = lm.split('lmp')[0]+'poscar'
    print(vsp)
    if 'conf.poscar' in os.listdir(md_i):
        os.remove(f'{md_i}/conf.poscar')
    os.system(f'ln -s {vsp} {md_i}/conf.poscar')
    # elems = os.popen(f'sed -n 6p {md_i}/conf.poscar').read().split()
    # print(elems)
#     if 'C' and 'H' and 'N' and 'O' in elems and len(elems)==4:
#         print('yes, I m in the loop')
    lmp_inp_file = lmp_12
    #elif 'Ti' and 'C' and 'H' in elems and len(elems)==3:
    #    lmp_inp_file = lmp_123
    os.system(f'echo {lm} > {md_i}/lmp_path.txt')
    md_i = f'{md_path}/md.{str(i).zfill(2)}'
    os.system(f'cp {lmp_inp_file} {md_i}/in.lammps')
    os.system(f'cp {ref_lammps}/job.slurm {md_i}/')
    os.system(f'cp {lm} {md_i}/conf.lmp')
    if 'traj' not in os.listdir(md_i):
        os.mkdir(f'{md_i}/traj')

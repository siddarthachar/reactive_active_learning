import os
from ase.io import read, write
import numpy as np
md_path_main = os.getcwd()
accepted = np.loadtxt(f'{md_path_main}/accepted.interm',dtype=str)
accepted_md = ['md.'+x.split('.')[0] for x in accepted]
#lmp_files = sorted([f'{md_path_main}/struct.packed/'+x for x in os.listdir(f'{md_path_main}/struct.packed') if x in f'{accepted}'])

lmp_files = sorted([f'{md_path_main}/struct.packed/'+x.split('.')[0]+'.lmp' for x in os.listdir(f'{md_path_main}/struct.packed') if x in f'{accepted}'])

dp_path = '../train/'
for i in range(4):
    os.system(f'ln -s {dp_path}/dp{i}/graph.pb {md_path_main}/graph.{str(i).zfill(3)}.pb')


ref_lammps = '/bgfs/kjohnson/pbs13/DeePMD/reactive_active_learning/methyle_03/sample-lammps'
lmp_12 = f'{ref_lammps}/in-123.lammps' # CHNO

for i, lm in enumerate(lmp_files):
    if f'{accepted_md[i]}' not in os.listdir(md_path_main):
        os.mkdir(f'{md_path_main}/{accepted_md[i]}')
    md_i = f'{md_path_main}/{accepted_md[i]}'
    vsp = lm.split('lmp')[0]+'poscar'
    if 'conf.poscar' in os.listdir(md_i):
        os.remove(f'{md_i}/conf.poscar')
    os.system(f'ln -s {vsp} {md_i}/conf.poscar')
    lmp_inp_file = lmp_12
    os.system(f'echo {lm} > {md_i}/lmp_path.txt')
    os.system(f'cp {lmp_inp_file} {md_i}/in.lammps')
    os.system(f'cp {ref_lammps}/job.slurm {md_i}/')
    os.system(f'cp {lm} {md_i}/conf.lmp')
    if 'traj' not in os.listdir(md_i):
        os.mkdir(f'{md_i}/traj')

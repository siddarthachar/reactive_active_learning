
from ase.io import read
from ase.neb import NEB
from ase.optimize import FIRE
import numpy as np
from deepmd.calculator import DP
import os
from pathlib import Path


here=int(os.getcwd().split('/')[-1].split('.')[-1]) # this is neb.xx
# structpath = '/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/reaction/test/'

curr=Path().absolute()
main_path=str(curr.parent.parent.parent)
structpath = f'{main_path}/init.struct/'

react_no = int(os.getcwd().split('/')[-2])
initial = read(f'{structpath}/ReactantOnSlab-{react_no}.poscar')
final = read(f'{structpath}/add_adsorbates_mol_only_reorder/ProdOnSurface-{str(int(here)).zfill(3)}.poscar')

pick = np.random.randint(0,4)
curr=Path().absolute()
main_path=str(curr.parent.parent.parent.parent)
#Change dp path to relative
dp_path_pick = f'{main_path}/train/dp{pick}/graph.pb'

images = [initial]
for i in range(5):
    image = initial.copy()
    image.calc = DP(model=f'{dp_path_pick}')
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
qn = FIRE(neb, trajectory='neb.traj')
qn.run(fmax=0.5, steps=500)

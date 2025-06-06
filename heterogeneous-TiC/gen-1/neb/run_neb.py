
from ase.io import read
from ase.neb import NEB
from ase.optimize import FIRE
import numpy as np
from deepmd.calculator import DP
import os

# change these

here=int(os.getcwd().split('/')[-1])
structpath = os.getcwd()

initial = read(f'{structpath}/rlx-c2h6-h-ads-0.poscar')
final = read(f'{structpath}/add_adsorbates_mol_only_reorder/ads-{here}.poscar')

pick = np.random.randint(0,4)
dp_path_pick = f'/bgfs/kjohnson/ska31/6AA/reactive_active_learning/TiC-methane-coupling/gen-1/train/dp{pick}/graph.pb'

images = [initial]
for i in range(5):
    image = initial.copy()
    image.calc = DP(model=f'{dp_path_pick}')
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
qn = FIRE(neb, trajectory='neb.traj')
qn.run(fmax=0.5)

from ase.io import read, write
from ase.visualize import view
import sys

trj = sys.argv[1] # path to lammpstrj
elems = sys.argv[2] # elems e.g. "Ti C"

elems = elems.split()
atoms = read(trj)
sym_index = atoms.get_atomic_numbers()
sym_dict = {i+1:sym for i, sym in enumerate(elems)} # {'Ti':1, "C":2}
new_chem_symb = [sym_dict[s] for s in sym_index]
atoms.set_chemical_symbols(new_chem_symb)
view(atoms)

import os
import sys

remove_force = int(sys.argv[1])

# Enter SMILES that you don't want in MD:
reject_smiles = ['C=N','O']

# all xyz files
xyz_files = sorted([x for x in os.listdir('struct.interm')])

def remove_forces_from_xyz(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    # Extract atomic positions
    atomic_positions = [line.split()[:4] for line in lines[2:]]
    with open(output_file, 'w') as f:
        f.write(str(num_atoms) + '\n')
        f.write(lines[1])  # Properties line
        for pos in atomic_positions:
            f.write(' '.join(pos) + '\n')

if remove_force:
	print('entered')
	if 'struct.interm.onlypos' not in os.listdir():
		os.mkdir('struct.interm.onlypos')
	for x in xyz_files:
		input_file = f'struct.interm/{x}'
		output_file = f'struct.interm.onlypos/{x}'
		remove_forces_from_xyz(input_file, output_file)
	print('Removed forces from all XYZ files - saved to struct.interm.onlypos')

# Printing out SMILES and comparing with reactant smiles
accepted = []
acc_file = open('accepted.interm','w+')
for x in xyz_files:
	interm_smiles = os.popen(f'python ~/bin/xyz2mol.py struct.interm.onlypos/{x} --use-huckel').read().split('\n')[0]
	interm_smiles_split = interm_smiles.split('.')
	alreadyseen = sum([0 if x not in reject_smiles else 1 for x in interm_smiles_split]) # if 0 then all mols in interm are new. 
	if alreadyseen == 0:
		accepted += [x]
		acc_file.write(f'{x}\n')
		print(f'Accepted: {x}: {interm_smiles}')
acc_file.close()




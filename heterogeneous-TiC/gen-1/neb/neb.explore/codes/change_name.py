import os

for i in range(42):
	os.system(f'mv {i} neb.{str(i).zfill(3)}')

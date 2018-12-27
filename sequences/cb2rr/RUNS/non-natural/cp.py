import sys, os

for i in range(63,74):
	#os.system('cp qsub.sh %d'%i)
	os.system('cp test_openmm.py gpu.sh %d'%i)

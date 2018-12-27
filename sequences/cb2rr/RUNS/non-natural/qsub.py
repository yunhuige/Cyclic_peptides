import sys, os
for i in range(63,73):
	os.chdir('%d'%i)
	os.system('qsub gpu.sh')
	os.chdir('../')

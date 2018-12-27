import sys, os
for i in range(63):
	os.chdir('%d'%i)
	#os.system('qsub qsub.sh')
	os.system('qsub gpu.sh')
	os.chdir('../')

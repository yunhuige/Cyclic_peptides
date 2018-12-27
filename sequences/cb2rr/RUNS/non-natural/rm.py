import sys, os
for i in range(63,74):
	os.chdir('%d/'%i)
	os.system('rm *trr')
	os.chdir('../')

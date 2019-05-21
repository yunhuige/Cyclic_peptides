import sys, os

for i in range(63,76):
    os.chdir('%d'%i)
    os.system('rm *#')
    os.system('rm *trr')
    os.system('rm *xtc')
    os.system('rm step*')
    os.system('rm *tpr')
    os.system('rm *edr')
    os.system('rm ligand.e*')
    os.chdir('../')

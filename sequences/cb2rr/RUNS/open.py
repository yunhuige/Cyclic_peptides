import sys, os

for i in range(63):
    os.chdir('%d'%i)
    os.system('open protein.gro')
    os.chdir('../')

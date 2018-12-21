import sys, os

for i in range(63):
    print i
    os.chdir('%d'%i)
    os.system('bash new_runme.sh')
    os.chdir('../')

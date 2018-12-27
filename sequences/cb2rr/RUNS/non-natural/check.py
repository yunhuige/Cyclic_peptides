import sys, os

for i in range(63,76):
    print i
    print os.path.isfile('%d/solvent_ions_equilibrated2.gro'%i)

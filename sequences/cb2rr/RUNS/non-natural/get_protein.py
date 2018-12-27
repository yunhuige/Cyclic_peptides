import sys, os
for i in range(63,74):
    os.chdir('%d'%i)
    os.system('gmx editconf -f solvent_ions_equilibrated2.gro -n index.ndx -o protein.gro')
    os.chdir('../')

import sys, os
for i in range(63):
    os.chdir('%d'%i)
    os.system('echo "1\n" | gmx editconf -f solvent_ions_equilibrated2.gro -n index.ndx -o protein.gro')
    os.chdir('../')

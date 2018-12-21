import sys, os

for i in range(63,73):
    print i
    os.chdir('%d'%i)

    with open('ligand_GMX.top','r') as file:
        data = file.readlines()
        print data[-6:]
        append = '; Include water topology\n#include "./amber99sbnmr1-ildn.ff/tip3p.itp"\n'
        append2 = '#ifdef POSRES_WATER\n; Position restraint for each water oxygen\n[ position_restraints ]\n;  i funct       fcx        fcy        fcz\n   1    1       1000       1000       1000\n#endif\n'
        append3 = '; Include topology for ions\n#include "./amber99sbnmr1-ildn.ff/ions.itp"\n'
        with open('new.top','w') as file:
            file.writelines(data[:-6])
            file.writelines('\n')
            file.writelines(append)
            file.writelines('\n')
            file.writelines(append2)
            file.writelines('\n')
            file.writelines(append3)
            file.writelines('\n')
            file.writelines(data[-6:])
    os.chdir('../')

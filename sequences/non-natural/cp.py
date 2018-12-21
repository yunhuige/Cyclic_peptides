import sys, os


for i in range(63,76):
    os.system('cp -r amber99sbnmr1-ildn.ff mdp spc216.gro %d'%i)

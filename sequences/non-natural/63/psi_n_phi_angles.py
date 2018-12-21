#!/usr/bin/env python

#usage='''

#chimera --start "Build Structure" --nogui --script "psi_n_phi_angles.py"\n

#**NOTE** You must hardcode your mutation info (look inside script)
#1) Loads pdb or gro into chimera
#2) Changes the psi and phi angles of residues that user specifies
#3) Saves new pdb

#'''
usage='''
This is a Chimera script that will...

Run in the terminal as:
chimera --nogui --script "psi_n_phi_angles.py ligand_GMX.gro"
'''

print usage
#chimera --script "script.py [pdb] [mutation] [residue.chain]"\n

# Import Modules:{{{
import os,sys,glob
import chimera
from Midas.midas_text import makeCommand as mc
import Midas as m
from chimera import replyobj
# }}}

# Methods:{{{

# Rotate Phi and Psi angles():{{{
files = glob.glob(str(sys.argv[1]))
#print files
#sys.exit()
sel=[':%d'%i for i in range(1,9)]
#print sel
#sys.exit()
psi = [111.376,136.591,113.452,-122.734,5.147,165.092,131.857,78.713,-133.349]
phi = [-95.475,-89.714,-151.946,50.611,-80.847,-131.562,-132.023,-120.535,58.513]
#    model = chimera.openModels.list()[int(model)]
mc('open %s'%files[0])
#sys.exit()
om = chimera.openModels
mlist = om.list()
#sys.exit()
model = str(mlist[0])
#sys.exit()
for i in range(1,9):
    print i
    mc('setattr r phi %s %s'%(str(phi[i-1]),sel[i-1]))
#    mc('setattr r psi %s %s'%(str(psi[i-2]),sel[i-2]))
mc('sel :')
mc('write format pdb selected %s %s'%(model,'test.pdb'))
print('write format pdb selected %s %s'%(model,'test.pdb'))
print '\n\n'
mc('close all')





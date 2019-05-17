#!/usr/bin/env python

usage='''

chimera --nogui --script "align.py [arg] [arg]"\n

**NOTE**
1) Loads pdb or gro into chimera
2)
3) Saves new pdb of mutated residues

'''
print(usage)

# Import Modules:{{{
import os,sys,glob
import chimera
import resCode
from phipsi import *
from chimera import runCommand as rc
from Midas.midas_text import makeCommand as mc
from chimera import replyobj
import Midas as m
import numpy as np
#import BuildStructure as BS
from chimera.baseDialog import ModalDialog
from chimera import openModels, Molecule
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError

# }}}

## Control:{{{
#def run_cmd(cmd, testing=False):
#    """Execute command-line command."""
#    print '>>', cmd
#    if not testing:
#        os.system(cmd)
#
#if len(sys.argv) >= 4:
#    run_cmd("echo '%s'"%usage)
#    sys.exit(1)
#
## Specific parameters:
#TESTING=False
#
## }}}
#
# Methods:{{{

def getChiAngles(r):
    return r.chi1,r.chi2,r.chi3,r.chi4

def getPhiAngles(r):
    return getPhi(r)

def getPsiAngles(r):
    return getPsi(r)



def rotate_Chi_angles(refID,mutID):
    """ This function will perform an optimization of rotamers with the
    great probability and select the rotamer that fits the Chi angles the
    best."""

    # Reference
    mc('sel %s'%refID)
    selRefRes = chimera.selection.currentAtoms()[0].residue
    print(selRefRes)
    refChiAngles = getChiAngles(selRefRes)
    print("The reference Chi angles are:\n %s"%str(refChiAngles))
    # Mutant
    mc('sel %s'%mutID)
    selMutRes = chimera.selection.currentAtoms()[0].residue
    print(selMutRes)
    mutChiAngles = getChiAngles(selMutRes)
    print("The mutant Chi angles are:\n %s"%str(mutChiAngles))
    #mc('setattr p chi %i sel'%(int(refChiAngles[0])))
    x = 1
    for i in range(len(refChiAngles)):
        if refChiAngles[i] != None:
            setChi(res=selMutRes,
                    chi=refChiAngles[i],
                    chiNum=x)
        x += 1
    mutChiAngles = getChiAngles(selMutRes)
    print("The mutant Chi angles are now:\n %s"%str(mutChiAngles))


def rotate_Phi_angles(refID,mutID):
    mc('sel %s'%refID)
    selRefRes = chimera.selection.currentAtoms()[0].residue
    print(selRefRes)
    mc('sel %s'%mutID)
    selMutRes = chimera.selection.currentAtoms()[0].residue
    print(selMutRes)
    refPhiAngles = getPhi(selRefRes)
    print("The reference Phi angles are:\n %s"%str(refPhiAngles))
    setPhi(res=selMutRes,phi=refPhiAngles)
    mutPhiAngles = getPhiAngles(selMutRes)
    print("The mutant Phi angles are now:\n %s"%str(mutPhiAngles))


def rotate_Psi_angles(refID,mutID):
    mc('sel %s'%refID)
    selRefRes = chimera.selection.currentAtoms()[0].residue
    refPsiAngles = getPsi(selRefRes)
    print("The reference Psi angles are:\n %s"%str(refPsiAngles))
    mc('sel %s'%mutID)
    selMutRes = chimera.selection.currentAtoms()[0].residue
    print(selMutRes)
    setPsi(res=selMutRes,psi=refPsiAngles)
    mutPsiAngles = getPsiAngles(selMutRes)
    print("The mutant Psi angles are now:\n %s"%str(mutPsiAngles))


# }}}

# Main:{{{


"""
1) center pdb - [crystal structure]
2) load - protein [crystal structure]
3) load - ligand [mutant]
4) use the "match" commmand to overlay [mutant] ontop the [crystal structure] ligand
5) lastly, remove the [crystal structure] ligand from the overlay structure
"""

ref = sys.argv[1]
print('\nLoading Crystal Structure...\n')
mc('open %s'%ref)

#TODO: automate the process of many mustant ligands (loop)
mutant = sys.argv[2]
print('\nLoading Mutant Ligand...\n')
mc('open %s'%mutant)
om = chimera.openModels
models = om.list()

print("###############################################################")
pdb_info = ['TITLE']#,'REMARK']
for m in openModels.list(modelTypes=[Molecule]):
    for header in pdb_info:
        try:
            for line in m.pdbHeaders[header]:
	            print('Model # = %s;\n Name = %s;\n%s;\n\n'%(
                    str(m.name),str(m.id),str(line)))
        except KeyError:
            print("Model # = %s;\n Name = %s;\n Has no header... \n\n"%(
                str(m.name),str(m.id)))
print("###############################################################\n\n")

print("There are %s models."%len(models))

# list of residues to match:
mod1 = ['#1:1','#1:3','#1:4']
mod0 = ['#0:21.B','#0:23.B','#0:24.B']

# Second Part: Match up the ligand with the crystal structure ligand:
#mc('mmaker %s %s'%(models[0],models[1]))
#mc('mmaker %s:.B %s pair ss alg smith show false iter 5.0'%(models[0],models[1]))
print("###############################################################")
mc('show %s:.B'%models[1])

# Part 2: Use Rotamers to alter the Chi angles to the specified

for i in range(0,len(mod0)):
    rotate_Chi_angles(refID=mod0[i],mutID=mod1[i])
    #rotate_Phi_angles(refID=mod0[i],mutID=mod1[i])
    #rotate_Psi_angles(refID=mod0[i],mutID=mod1[i])

mc('~sel')

# Second Part: Match up the ligand with the crystal structure ligand:
mc('mmaker %s:.B %s pair ss alg smith show false iter 15.0'%(models[0],models[1]))
#mc('mmaker %s %s'%(models[1],models[0]))

mc('sel #0:.B')
mc('show sel')


        ## Combine Models and Close Old Models
        #mc('combine #1,#0 modelID 0 close True')
        # Change the Chain ID
        #From = 'het'
        #To = SN[-1]
        #mc('changechains %s %s'%(From, To))
        #mc('delete %s'%spec1)
        # Save the new pdbs
        #outFile = pdb.replace('.pdb','_%s.pdb'%seqNames[x])
        #mc('write format pdb selected %s %s'%(model,outFile))
        #outFile = pdb.replace('.pdb','_%s.mol2'%seqNames[x])
        #mc('write format mol2 selected %s %s'%(model,outFile))
        #print('write format pdb selected %s %s'%(model,outFile))
        #print('write format mol2 selected %s %s'%(model,outFile))
        #print '\n\n'
        #mc('close all')


# }}}





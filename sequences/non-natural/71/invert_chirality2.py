#!/usr/bin/env python

usage='''
This is a Chimera script that will...

Run in the terminal as:
chimera --nogui --script "script_name.py"
'''
print(usage)

# Importing Modules:{{{
import os,sys,glob
import chimera
from Midas.midas_text import makeCommand as mc
import Midas as m
from chimera import replyobj
# }}}

# Selection
# What is the selection (Atom specific) for invert
# example:   Residue.Chain@Atom
#            17.B@CB
sel = ':3@CA'

# Main:{{{

files = glob.glob(str(sys.argv[1]))
for f in files:
    mc('open %s'%f)
    om = chimera.openModels
    mlist = om.list()
    model = str(mlist[0])

    # Now, calling on the invert function:
    mc('invert %s'%sel)

    if f.endswith('.pdb'):
        outFile = f.replace('.pdb','_invert2.pdb')
    elif f.endswith('.gro'):
        outFile = f.replace('.gro','_invert2.pdb')
    else:
        print('This file does not end with .pdb or .gro')
        exit(1)

    mc('sel :')
    mc('write format pdb selected %s %s'%(model,outFile))
    print('write format pdb selected %s %s'%(model,outFile))
    print '\n\n'
    mc('close all')


# }}}






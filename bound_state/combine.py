#!/usr/bin/env python

usage="""

chimera --nogui --script "combine.py"\n
"""

print(usage)
# Import Modules:{{{
import sys,os
import chimera
from Midas.midas_text import makeCommand as mc
from chimera import openModels, Molecule
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError

#if (sys.argv[1:] == ("-h" or "--h" or "--help" or "help")) or (len(sys.argv) >= 4):
#    print(usage)
#    exit(1)
# }}}



def load_pdb(pdb, verbose=False):
    """Loads a structure file into chimera."""

    print('\n\n\nLoading Crystal Structure...\n')
    mc('open %s'%ref)

    om = chimera.openModels
    models = om.list()

    if verbose:
        # Get the information from the pdb file loaded:
        print("###############################################################")
        pdb_info = ['TITLE']#,'REMARK']
        for m in openModels.list(modelTypes=[Molecule]):
            for header in pdb_info:
                try:
                    for line in m.pdbHeaders[header]:
                        print('Model # = %s;\n Name = %s;\n%s;\n\n'%(
                            str(m.id), str(m.name), str(line)))
                except KeyError:
                    print("Model # = %s;\n Name = %s;\n Has no header... \n\n"%(
                        str(m.id), str(m.name)))
        print("There are %s models."%len(models))
        print("###############################################################\n\n")
    return models



# MAIN:

ref = "Path to 1st pdb" # sys.argv[1]
ref1 = "Path to 1st pdb" # sys.argv[2]



models = load_pdb(ref) # model #1
models = load_pdb(ref1) # model #2
#mc('write format pdb %s %s'%(models[0],outFile))

# will save all models
mc('write format pdb %s %s'%(models,outFile))
mc("close all")






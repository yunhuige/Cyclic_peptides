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
    mc('open %s'%pdb)

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

models = load_pdb(sys.argv[1]) # model #1
mc('open %s'%sys.argv[2])
mc("combine #0,1")
outFile = "new.pdb"
om = chimera.openModels
models = om.list()
# will save all models
mc('write format pdb %s %s'%(models[2], outFile))
mc("close all")






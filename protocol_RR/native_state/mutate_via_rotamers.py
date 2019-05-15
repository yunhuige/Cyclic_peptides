#!/usr/bin/env python

usage='''

chimera --nogui --script "mutate_via_rotamers.py [arg] [arg]"\n

This program will use the rotamer library to alter the given structure.
The user has the option to mutate the residues or just adjust the torsions.
Required Inputs:
    - Reference Structure (pdb or gro)
Optional Inputs:
    - List of sequences that will result in the mutated structures
Output:
    - pdb or gro file
'''
# Import Modules:{{{
import sys
import chimera
from Midas.midas_text import makeCommand as mc
from chimera import openModels, Molecule
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError

if (sys.argv[1:] == ("-h" or "--h" or "--help" or "help")) or (len(sys.argv) >= 4):
    print(usage)
    exit(1)
# }}}

# Methods:{{{

def mutate_from_rotamers(refID, resType=None, replace=False, verbose=False):
    """ This function will perform an optimization of rotamers with the
    great probability and select the rotamer that fits the Chi angles the
    best.

    :param str refID: the chimera residue ID for the reference
    :param str resType: the residue that will result after mutation
    """

    # Residues are uppercase
    AA = {'R' : 'ARG','H' : 'HIS','K' : 'LYS','D' : 'ASP',
            'E' : 'GLU','S' : 'SER','T' : 'THR','N' : 'ASN',
            'Q' : 'GLN','C' : 'CYS','U' : 'SEC','G' : 'GLY',
            'P' : 'PRO','A' : 'ALA','V' : 'VAL','I' : 'ILE',
            'L' : 'LEU','M' : 'MET','F' : 'PHE','Y' : 'TYR',
            'W' : 'TRP'}

    execute = True
    if resType != None:
        # Make sure that the string is uppercase:
        resType = resType.upper()
        # If given a 1 letter symbol, then get the 3 letter identity:
        if len(resType) != 3:
            resType = AA[str(resType)]

    # Reference
    mc('sel %s'%refID)
    selR = chimera.selection.currentAtoms()[0].residue
    selRType = selR.type
    if resType == selRType:
        if replace == False:
            print("pass\n\n")
            execute = False

    if execute:
        # Borrowing (The code was rewritten) the next 15 lines of code from:
        # Computational and Visualization Techniques for Structural Bioinformatics ...
        # By Forbes J. Burkowski
        # Page 269
        rotIx = 0
        try:
            (flag, rots_L) = getRotamers(res=selR, resType=resType)
            numRots = len(rots_L)
            #print("\n")
            for i in range(numRots):
                resRot = rots_L[rotIx].residues[0]
                if verbose:
                    print("%s %s\n%s %s"%("Rotamer with index: ", rotIx, \
                          "probability ", rots_L[rotIx].rotamerProb))
                if rotIx == 0:
                    print("Using rotamer with highest probability...")
                    print("%s %s\n%s %s"%("Using rotamer with index: ", rotIx, \
                          "probability ", rots_L[rotIx].rotamerProb))
                    useRotamer(selR, [rots_L[rotIx]])
                rotIx += 1
        except NoResidueRotamersError:
            print(selR.id.position, selR.type, "No rotamers")

# }}}

# Main:{{{
# Reference
ref = sys.argv[1]
print('\nLoading Crystal Structure...\n')
mc('open %s'%ref)

om = chimera.openModels
models = om.list()

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

# TODO: Get rid of these hardcodes if possible. Otherwise they need to be moved
# and mentioned in the usage message.
# NOTE: List of residues that should not be mutated if they are in their correct
# positions according to the sequence.
ref_seq = ['#0:21.B','#0:22.B','#0:23.B','#0:24.B','#0:25.B','#0:26.B',
        '#0:27.B','#0:28.B']
ref_res = ["PHE","GLU","6CW","LEU","ASP","TRP","GLU","PHE"]


# Open the file containing the sequences and number of each of the sequences
seq_number, new_seq = [], []
if sys.argv[2]:
    with open(str(sys.argv[2]), "r") as file:
        for line in file.readlines():
            L = [i.strip("\n") for i in str(line).split(" ")]
            seq_number.append(L[0])
            new_seq.append(L[1])

# Part 2: Use Rotamers library:
nth_seq = 0
for seq in new_seq:
    nth_res = 0
    for new_res in seq:
        mutate_from_rotamers(refID=ref_seq[nth_res],resType=new_res)
        nth_res += 1

    outFile = "2axi_seq_%s.pdb"%(seq_number[nth_seq])
    print('write format pdb %s %s'%(models[0],outFile))
    mc('write format pdb %s %s'%(models[0],outFile))
    # Stop here to check
    # comment or uncomment the break
    break
    nth_seq += 1

mc('~sel')
mc('sel #0:.B')
mc('show sel')

# }}}





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
import sys,os
import chimera
from Midas.midas_text import makeCommand as mc
from chimera import openModels, Molecule
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError

if (sys.argv[1:] == ("-h" or "--h" or "--help" or "help")) or (len(sys.argv) >= 4):
    print(usage)
    exit(1)
# }}}

# Methods:{{{
def run_cmd(cmd, testing=False):
    """Execute command-line command."""

    if testing:
        print('>>', cmd)
    else:
        os.system(cmd)

def replace_text_in_file(FILE, text, replacement, out):
    with open(FILE, "r") as file:
        with open(out, "w") as f:
            f.seek(0)
            for line in file.readlines():
                if text in line:
                    f.write(line.replace(text,replacement))
                else:
                    f.write(line)
            f.truncate()
            f.close()
        file.close()


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


def swapAA(resName,ID):
    '''Mutates an amino acid to a new specified residue.'''

    # Mutating amino acids:
    print('Mutation: ',ID,' to ',resName)
    # Swap a residue of choice for another  amino acid:
    mc('swapaa %s %s'%(resName,ID))


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
            print("Chimera selection: %s, %s pass\n\n"%(selR,selRType))
            execute = False

    if execute:
        print("\nChimera selection: %s, %s --> %s"%(selR,selRType,resType))
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
            #NOTE: Force the mutation #TODO: fix this with a conditional statement
            print("Forcing a mutation...")
            print("mutating %s --> %s"%(selR.type,resType))
            swapAA(resType,refID)

# }}}

# Main:{{{
# Reference
ref = sys.argv[1]

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
for i in range(len(new_seq)):
    # Open the pdb in chimera
    models = load_pdb(ref)

    print("##########################################")
    print("%s, %s"%(seq_number[i],new_seq[i]))
    print("##########################################")
    for j in range(len(new_seq[i])):
        mutate_from_rotamers(refID=ref_seq[j],resType=new_seq[i][j])
    outFile = "2axi_seq_%s.pdb"%(seq_number[i])
    # Write the first pdb file with protein and ligand
    print('write format pdb %s %s'%(models[0],outFile))
    mc('write format pdb %s %s'%(models[0],outFile))



    #NOTE: This is only to change DRP --> PRO
    # If the file contains residues named "DPR", then we need to change them to "PRO"
    newoutFile = "2axi_.pdb"
    replace_text_in_file(FILE=outFile, text="DPR", replacement="PRO", out=newoutFile)
    #replace_text_in_file(FILE=outFile, text="HIS", replacement="HIE", out=newoutFile)

    #NOTE:
    #NOTE: Everything included in the ####### after this point is added
    # with the intent to separate the ligand and protein.
    ## Here, we will separate and create a new pdb of only the ligand, then we
    # will add hydrogens.
    # Write the 2nd and 3rd pdbs which are  protein and ligand separated
###############################################################################

    mc("addh")
    mc("sel :.B")
    #FIXME
    lig_outFile = "2axi_ligand.pdb"
    print('write format pdb selected %s %s'%(models[0],lig_outFile))
    mc('write format pdb selected %s %s'%(models[0],lig_outFile))

    mc("sel :.A")
    pro_outFile = "2axi_protein.pdb"
    print('write format pdb selected %s %s'%(models[0],pro_outFile))
    mc('write format pdb selected %s %s'%(models[0],pro_outFile))

    out1 = "ligand.pdb"
    replace_text_in_file(FILE=lig_outFile, text="DPR", replacement="PRO", out=out1)
    #replace_text_in_file(FILE=outFile, text="HIS", replacement="HIE", out=newoutFile)
    out2 = "protein.pdb"
    replace_text_in_file(FILE=pro_outFile, text="DPR", replacement="PRO", out=out2)
    #replace_text_in_file(FILE=outFile, text="HIS", replacement="HIE", out=newoutFile)

    if not os.path.exists(str(i)):
        run_cmd("mkdir %s"%str(i)) # make a directory
    run_cmd("mv %s %s"%(outFile,i))
    run_cmd("mv %s %s"%(out1,i))
    run_cmd("mv %s %s"%(out2,i))
    run_cmd("mv %s %s"%(newoutFile,i))
    run_cmd("mv %s %s"%(lig_outFile,i))
    run_cmd("mv %s %s"%(pro_outFile,i))


###############################################################################


    mc('close all')
    # Close the pdb after it has been written and go to the next iteration.
    # Stop here to check
    # comment or uncomment the break
    #break



# }}}

#exit(1)
## Open up the mapping.txt file:
#
#mapping = "/Users/tuc41004/github/Cyclic_peptides/protocol_RR/old_mapping.txt"
#
#with open(mapping, "r") as file:
#    for line in file.readlines():
#        From = line.split(" --> ")[0]
#        To = line.split(" --> ")[1]
#        # NOTE: In this case we are mapping in reverse (To --> From)
#        if not os.path.exists(From):
#            run_cmd("mkdir %s"%From) # make a directory
#
#        FILES = ["2axi_seq_%s.pdb"%(To), "2axi_seq_%s_ligand.pdb"%(To), "2axi_seq_%s_protien.pdb"%(To)]
#
#        for i in range(len(FILES)):
#            new_dir = "./"+From+"/"+FILES[i].replace(To,From)
#            run_cmd("mv %s %s"%(FILES[i],new_dir))
#        #run_cmd("")
#        #run_cmd("")
#








#!/user/bin/env python

# Python Libraries:{{{
import numpy as np
import sys, os, platform
#}}}

# Control/Functions:{{{
TESTING = False #True 
def run_cmd(cmd, testing=TESTING):
    """Execute command-line command."""
    if not testing:
        os.system(cmd)
    else:
        print '>>', cmd

def edit_top(top_in, top_out):
        """Given a gmx top file, add necessary include lines and delete defaults section"""

        fin = open(top_in, 'r')   # open th3 file handle
        # read file as a list of strings
        lines = fin.readlines()

        fin.close()   # close th4 file handle

        # Let's look at the first 10 lines
        for i in range(10):
            print 'Line %d:'%i, lines[i].strip('\n')

        # let's get rid of lines 2,3,4
        for i in range(3):
            lines.pop(2)

        # Next, put an include statement in line 1 
        lines.insert(1,'# include "amber14sb.ff/forcefield.itp"\n') 

        # put two include statements in at line -9
        lines.insert(-6,'# include "amber14sb.ff/tip3p.itp"\n')
        lines.insert(-6,'# include "amber14sb.ff/ions.itp"\n')

        ## Write the modified top file 
        fout = open(top_out,'w')
        fout.writelines(lines)                                                    

#}}}

# Hardcode: User Imputs:{{{
# Using tleap to build the structure and top file and using acpype to convert
# them to GROMACS format

# Directory for Chimera app
chimera='/Applications/Chimera.app/Contents/MacOS/chimera'
# Directory for acpype scipt
acpype='/Users/allisonterry/voelzlab/custom_topologies/scripts/acpype.py'
# tleap *.in files
in_files='/Users/allisonterry/Cyclic_peptides/allison/tleap_inputs/'
# Provide the working directory for this program
WD='/Users/allisonterry/Cyclic_peptides/allison/'
# Provide the output RUNS Directory
RUNS=WD+'RUNS/'

# Provide a list of the sequences
seq = [1,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
        25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
        48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
        67,68,69,70,71,72,73,74,75,76,77,78,79,80]


# List the files that need to be moved over to the RUNS directories
neededFiles_dir = '/Users/allisonterry/Cyclic_peptides/allison/needed_files/'
neededFiles = ['amber14sb.ff','mdp']

# Create a file for mapping information
logFile = "mapping.txt"
if os.path.isfile(logFile):
    run_cmd('mv mapping.txt old_mapping.txt')
run_cmd("touch %s"%logFile)

# Create an error file:
errorFile = 'errors.txt'
if os.path.isfile(errorFile):
    run_cmd('mv errors.txt old_errors.txt')
run_cmd("touch %s"%errorFile)

# Create a RUNS Output directory if one doesn't exist
if not os.path.isdir('%s'%(RUNS)):
    run_cmd('mkdir %s'%(RUNS))

# Clean up these files:
acpypeFiles = ['ligand.prmtop','ligand.crd']
for f in acpypeFiles:
    run_cmd('rm -v %s%s'%(WD,f))


#}}}

# Main:{{{
# Loop through each of the tleap *.in files and build the structure files
# Prepare the RUN directories for the simulations
n=0
for i in range(0,len(seq)):
    # Make the Directory if it doesn't exist
    if not os.path.isdir('%s%s'%(RUNS,i)):
        run_cmd('mkdir %s%s'%(RUNS,i))

    # Figure out if we have this sequence...
    if os.path.isfile('%stleap_%s.in'%(in_files,seq[i])):
        print 'Reading file: tleap_%s.in ...'%(seq[i])

        # Use tleap to read a specific tleap_*.in file for the ith ligand
        run_cmd('tleap -f %s/tleap_%s.in'%(in_files,seq[i]))

        # Use acpype to generate gromacs structure and topology files
        run_cmd('%s -p ligand.prmtop -x ligand.crd -b ligand;'%acpype)
        # Name the files produced from acpype:
        gro = 'ligand_GMX.gro'
        top = 'ligand_GMX.top'

        # Creating variable for output directory (.../RUNS/#)
        outDir = str(RUNS)+str(n)

        # Copy the important files over to its corresponding directory.
        for j in range(len(neededFiles)):
            run_cmd('cp -rv %s%s %s'%(neededFiles_dir,neededFiles[j],outDir))
        run_cmd('mv -v %s %s'%(gro,outDir))
        run_cmd('mv -v %s %s'%(top,outDir))

        for f in acpypeFiles:
            run_cmd('mv -v %s%s %s'%(WD,f,outDir))

        # Edit the topology file
        top = outDir+'/'+str(top)
	out = outDir+'/top_edit.top'
        edit_top(top, out)

        # Using Chimera to modify the chirality of a particular ligand
        run_cmd('%s --nogui --script invert_chirality.py %s/ligand_GMX.gro'%(chimera,outDir))

        # Using Chimera to modify dihedral angles of a particular ligand
        run_cmd('%s --nogui --script psi_n_phi_angles.py %s/ligand_GMX_invert.pdb'%(chimera,outDir))

        # Check to make sure that psi_n_phi_angles.py script worked correctly:
        if not os.path.isfile('%s/ligand_GMX_invert_angles.pdb'%(outDir)):
            print("NOTE: Missing File. Appending to %s..."%errorFile)
            fout = open(errorFile,"a+")
            fout.write("%s\n"%(seq[i]))
            fout.close()

        #NOTE: Remove the old files and overwrite old name to be consistent
        run_cmd('cp -v %s %s'%(outDir+"/ligand_GMX_invert_angles.pdb",
            outDir+"/ligand_GMX_invert.pdb"))
        run_cmd('rm %s'%(outDir+'/ligand_GMX_invert_angles.pdb'))

        # Copy the runme.sh for simulation script to each RUN
        run_cmd('cp runme.sh %s'%(outDir))

    # Append to a file that contains the mapping from RUN# --> sequence#
    print "Appending to %s..."%logFile
    fout = open(logFile,"a+")
    fout.write("%s --> %s\n"%(n,seq[i]))
    fout.close()
    n += 1


#}}}




#!/usr/bin/env python

import os, sys, glob, commands, argparse, tempfile
from Mol2File import *

class Parameters():
    pass #IE do nothing

# Functions

def run_cmd(cmd, Verbose=True):
    #Executes a command-line cmd and prints the output, unless specified not to.
    #Could use call for this rather than getoutput
    print '>>', cmd
    output = commands.getoutput(cmd)
    if Verbose: print output
    return


def build_antechamber(InputFileName, ResidueName, CappingAtoms, CNTerm, Forcefield, NetCharge,SmearCharge, HeadAtom='N', TailAtom='C'):

   
    print 'Building from file', InputFileName, '...'

    # Check to see if it's a mol2 or pdb file
    if InputFileName[-4:] == '.pdb':
        InputFileType = 'pdb'
        in_mol2  = InputFileName.replace('.pdb','.preedit.mol2') # Will be InputFileName.preedit_mol2
        out_mol2 = InputFileName.replace('.pdb','.edited.mol2') # Will be InputFileName.edited.mol2
    elif InputFileName[-5:] == '.mol2':
        InputFileType = 'mol2'
        in_mol2  = InputFileName.replace('.mol2','.preedit.mol2') # Will be InputFileName.preedit_mol2
        out_mol2 = InputFileName.replace('.mol2','.edited.mol2') # Will be InputFileName..edited.mol2
    else:
        print 'Input file must be either a *.pdb or *.mol2 file!  Exiting...'
        sys.exit(1)

    # make mol2 for residue <resname>

    if (1):  # we can skip this step if we just want to re-run  the charge modifications

        if (Forcefield=='Gaff') or (Forcefield=='G') or (Forcefield=='gaff'):
            # run antechanber to get a mol2 file with GAFF atomtypes and  BCC charges
            #cmd =  'antechamber -i %s -fi mol2 -o %s -fo mol2 -at gaff -nc %d -j 5'%(mol2file, in_mol2, NetCharge) 
            #print '>>', cmd
            #run_cmd( cmd )"""
            command='antechamber -i %s -fi %s -o %s -fo mol2 -at gaff -nc %d -j 5 -c bcc'%(InputFileName, InputFileType, in_mol2, NetCharge)
            run_cmd(command)
            #call([command])

        else:  # use AMBER
            # run antechanber to get a mol2 file with AMBER atomtypes and  BCC charges
            #run_cmd( 'antechamber -i %s -fi pdb -o %s -fo mol2 -at amber -nc %d -j 5'%(mol2file, in_mol2, NetCharge) )
            command='antechamber -i %s -fi %s -o %s -fo mol2 -at amber -nc %d -j 5 -c bcc'%(InputFileName, InputFileType, in_mol2, NetCharge)
            run_cmd(command)
            #jcall([command])
            


    # We'll keep the charges produced by GAFF (no overwriting with the AMBER partial charges),
    # but for the capping atoms that will get removed, set the charge to zero (these will later get removed when building the resideu in tleap)

    ResObject=Mol2File(in_mol2)#Casts the mol2file object type to the residue so that the edit_charge and net_charge methods can be called on it

    # 1) set all the capped atoms to zero charge 
    for EveryAtom in CappingAtoms:
       # m.edit_charge(EveryAtom, 0.0)
       ResObject.edit_charge(EveryAtom, 0.0)
    print 'net charge after zeroing capped residues:%f' % (ResObject.net_charge())  #SHould this be zero at this point?  Seems to be at .0022 after zeroing capped residues, seemingly negligible

    # 2_ Then, shift all the atomic charges to get an integer net charge
    #m.make_net_charge_integer(NetCharge, exclude_names=capping_atoms)
    #print 'net charge after shifting to net charge:', m.net_charge()
    
    ResObject.make_net_charge_integer(SmearCharge, exclude_names=CappingAtoms)
    print 'net charge after shifting to net charge:%f' % (ResObject.net_charge())

    print 'Writing to', out_mol2
    ResObject.write(out_mol2)

    print ResObject

    return out_mol2


def tleap_residue_lines(res_mol2file, res, CNTerm, capping_atoms, HeadAtom='N', TailAtom='C'):
    #Returns a string corresponding to the tleap lines needed to build a residue (i.e. remove the terminal caps, and specify connect points)

    if not os.path.exists(res_mol2file):
        print '***WARNING!*** Cannot find residue mol2file "%s"... Skipping.'%res_mol2file 
        return ''

    ResObject = Mol2File(res_mol2file)    

    leaptxt =  '%s = loadMol2 %s\n'%(res, res_mol2file)
    # Now built the terminal residue stuff using a similar scheme, with GAFF, etc as described in README; go ahead and try adding here as well.
    # %s = loadMol2 %s.edited.mol2\n%(res,res) Think this line is just a formatting issue

    if CNTerm == 'B':
        capping_indices = []  # ignore capping indices
    else:
        capping_indices = [ ResObject.get_atom_id_by_name(atomname) for atomname in capping_atoms] 

    leaptxt += "\n# Delete N-terminal acetyl and C-terminal dimethyl amine group\n"
    for i in range(len(capping_indices)): 
    	leaptxt += "remove %s %s.1.%d  # %s\n"%(res,res,capping_indices[i],capping_atoms[i])
    leaptxt += '\n'
  
    leaptxt += '# Should still be net neutral -- all removed atoms shoulve have had charge 0.0\n\n'

    leaptxt += 'set %s name "%s"\n'%(res,res)
    leaptxt += 'set %s.1 name "%s"\n'%(res,res)
    leaptxt += '\n'

    if (CNTerm == 'C') or (CNTerm == 'X'):
        leaptxt += '# Define head atom\n'
        leaptxt += 'set %s head %s.1.%d\n'%(res,res,ResObject.get_atom_id_by_name(HeadAtom))
        leaptxt += 'set %s.1 connect0 %s.1.%d\n'%(res,res,ResObject.get_atom_id_by_name(HeadAtom))
        leaptxt += '\n'

    if (CNTerm == 'N') or (CNTerm == 'X'):
        leaptxt += '# Define tail atom\n'
        leaptxt += 'set %s tail %s.1.%d\n'%(res,res,ResObject.get_atom_id_by_name(TailAtom))
        leaptxt += 'set %s.1 connect1 %s.1.%d\n'%(res,res,ResObject.get_atom_id_by_name(TailAtom))
        leaptxt += '\n'

    leaptxt += 'savepdb %s %s.pdb\n\n\n'%(res,res)
 
    return leaptxt





if __name__ == "__main__":     
	
#Step1 Code-----------------------------------------------------------------------
    
	import argparse

        p = Parameters()#In order for the argument variables to be used later in the script, they must be members of a class (p/parameters in this case)    
	parser = argparse.ArgumentParser(description='Arguments needed to build residue topologies')

        parser.add_argument('InputFileName', action="store",help='Name of the input file.  Can either be *.pdb or *.mol2')
	parser.add_argument('ResidueName', action="store",help='Name of the output <ResidueName>.off residue')
	parser.add_argument('CappingAtomsFile', action="store", help='Filename with the names of each capping atom, one per line')
	parser.add_argument('CNTerm',action="store",choices='NCXB', help="""Specify N-term (N), C-term (C), neither (X) or both (B, i.e. a complete moleculr that is not part of chain) Details:
# (N) an N-terminal residue; will only specify "tail" atom to connect; capping atoms are the C-term cap
# (C) a C-terminal residue; will only specify "head" atom to connect; capping atoms are the N-term cap
# (X) a normal residue; will specify both ""head" and "tail" atoms to connect, capping atoms both caps
# (B) a self-contained residue; no ""head" or "tail" atoms, there should be no capping atoms to remove
""")
	parser.add_argument('Forcefield',action="store",help='Specify the GAFF forcefield (G) or the Amber forcefield (A)')
	parser.add_argument('NetCharge',action="store",type=float,help='Specify the net charge of the residue')
	parser.add_argument('SmearCharge',action="store",type=float,help='Specify the charge need to smear for the whole system')
	#Residue, CappingAtoms, CNTerm, Forcefield, NetCharge

        parser.add_argument('-s', dest='HeadAtom', help='Name of head atom for connecting residues',
        default='N')
        parser.add_argument('-t', dest='TailAtom', help='Name of tail atom for connecting residues',
        default='C')

        #parser.add_argument('--HeadAtom', action='store_const', const='N', default='N')
        #parser.add_argument('--TailAtom', action='store_const', const='C', default='C')

	args =  parser.parse_args(namespace=p)

        #>>> parser = argparse.ArgumentParser()
        #>>> parser.add_argument('--foo')
        #>>> parser.parse_args(args=['--foo', 'BAR'], namespace=c)
        #>>> c.foo
        #'BAR'

        # read in cappinhg atoms from file
        if p.CNTerm == 'B':
            CappingAtoms = [] # ignore capping atoms file if we are building a complete molecule with no caps
        else:
            fin = open(p.CappingAtomsFile, 'r')
            CappingAtoms = fin.read().strip().split()
            fin.close()

	out_mol2 = build_antechamber(p.InputFileName, p.ResidueName, CappingAtoms, p.CNTerm, p.Forcefield, p.NetCharge, p.SmearCharge,HeadAtom=p.HeadAtom, TailAtom=p.TailAtom)
    

    #mol2files = glob.glob('residues/*.mol2')
    #mol2files.sort()
    #print 'mol2files', mol2files

    #for mol2file in mol2files:
        #build_antechamber(mol2file, UseGAFFAtomTypes=True)


#Step2 Code--------------------------------------------------------------------------


	# Use paramchk to generate a file with modified forcefield parameters we may need
	out_mol2_frcmod = out_mol2.replace('.mol2','.mol2.frcmod') 
	cmd = 'parmchk -i %s -f mol2 -o %s -p $AMBERHOME/dat/leap/parm/gaff.dat'%(out_mol2, out_mol2_frcmod)
	print '>>', cmd
	os.system(cmd)


#Step 3 Code--------------------------------------------------------------


        # Write the dlm_edit.frcmod file for extra missing parms
        dlm_text = """remark goes here - Originally created by DLM, VAV edited to mix-n-match AMBER (uppercase) and GAFF (lowercase) atomtypes
MASS

BOND
C -NT   486.3    1.340       DLM: GAff c2-n3
c -N    478.2    1.345       VAV: GAFF c -n  SOURCE1  1235    0.0162    0.0215
C -n    478.2    1.345       VAV: GAFF c -n  SOURCE1  1235    0.0162    0.0215
N -c    478.2    1.345       VAV: GAFF c -n  SOURCE1  1235    0.0162    0.0215
n -C    478.2    1.345       VAV: GAFF c -n  SOURCE1  1235    0.0162    0.0215
CT -ss  227.0    1.810	     GYH: copied from parm99.dat CT -S 

ANGLE
O -C-NT    74.1       122.26     DLM GAFF n3-c -o 
C -NT-CT    50.0      109.50     DLM analogy to parm99
C -NT-H    49.1      119.38      DLM GAFF c2-n3-hn
C -CT-NT    80.0      111.20     DLM analogy to parm99
CT-C -NT    64.6      112.26     DLM from gaff c3-c -n4; for whatever reason there is no parameter in gaff or parm99 for CT-C-NT or similar
O -C -n     75.8      122.03     VAV: GAFF n -c -o  SOURCE3  221   1.7197   2.3565
o -c -N     75.8      122.03     VAV: GAFF n -c -o  SOURCE3  221   1.7197   2.3565
C -n -c3    63.9      121.35     VAV: GAFF c -n -c3 SOURCE3   54   1.7456   2.3808
c -N -CT    63.9      121.35     VAV: GAFF c -n -c3 SOURCE3   54   1.7456   2.3808
c3-c -N     67.9      115.15     VAV: GAFF c3-c -n  SOURCE3  153   1.9677   2.7443
CT-C -n     67.9      115.15     VAV: GAFF c3-c -n  SOURCE3  153   1.9677   2.7443
c -N -H     49.2      118.46     VAV: GAFF c -n -hn SOURCE3  157   1.8119   2.4094
C -n -hn    49.2      118.46     VAV: GAFF c -n -hn SOURCE3  157   1.8119   2.4094
ca-c -N     63.0      120.00     VAV: coped from parm99.dat CA-C CA:  changed from 85.0  bsd on C6H6 nmodes; AA
H1-CT-ss    50.0      109.50	 GYH: copied from parm99.dat H1-CT-S
CT-CT-ss    50.0      114.70	 GYH: copied from parm99.dat CT-CT-S
CT-ss-c3    62.0       98.90	 GYH: copied from parm99.dat CT-S-CT
CT-ss-ss    68.0      103.70     GYH: copied from parm99.dat CT-S-S

DIHE
O -C-NT-CT     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
O -C-NT-H     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
C-NT-CT-CT     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
C-NT-CT-C     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
C-NT-CT-H1     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
CT-C-NT-CT     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
CT-C-NT-H     6    1.80          0.0             3.        Analogy to parm99 X-CT-NT-X
CA-CA-CA-CA   1    3.625       180.000           2.000      same as X -ca-ca-X 
CA-CA-CA-HA   1    6.650       180.000           2.000      same as X -c2-c2-X 
CA-CA-CA-CT   1    6.650       180.000           2.000      same as X -c2-c2-X 
CA-CA-CT-CT   1    0.000         0.000           2.000      same as X -c2-c3-X 
CA-CA-CT-HC   1    0.000         0.000           2.000      same as X -c2-c3-X 
CA-CT-CT-CT   1    0.156         0.000           3.000      same as X -c3-c3-X 
CA-CT-CT-HC   1    0.156         0.000           3.000      same as X -c3-c3-X 
C -CT-CT-CT   1    0.156         0.000           3.000      same as X -c3-c3-X 
C -CT-CT-HC   1    0.156         0.000           3.000      same as X -c3-c3-X 
CT-C -CT-CT   1    0.000       180.000           2.000      same as X -c -c3-X 
CT-C -CT-HC   1    0.000       180.000           2.000      same as X -c -c3-X 
CT-CA-CA-HA   1    6.650       180.000           2.000      same as X -c2-c2-X 
CT-CT-CT-HC   1    0.156         0.000           3.000      same as X -c3-c3-X 
CT-CT-C -O    1    0.000       180.000           2.000      same as X -c -c3-X 
O -C -CT-HC   1    0.000       180.000           2.000      same as X -c -c3-X 
HA-CA-CA-HA   1    6.650       180.000           2.000      same as X -c2-c2-X 
HC-CT-CT-HC   1    0.156         0.000           3.000      same as X -c3-c3-X 
O-C-n-c3      4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
CT-C-n-c3     4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
o-c-N-H       4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
O-C-n-hn      4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
c3-c-N-CT     4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
o-c-N-CT      4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
c3-c-N-H      4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
CT-C-n-hn     4   10.000       180.000           2.000      VAV analogy to gaff X -c -n -X:  AA,NMA (no c-n3, c-n4, c-nh)
ca-c-N-H      4   10.00        180.0             2.         VAV: copied from  parm99.dat X -C -N -X :  AA,NMA
ca-c-N-CT     4   10.00        180.0             2.         VAV: copied from  parm99.dat X -C -N -X :  AA,NMA
H1-CT-ss-c3   3    1.00          0.0             3.	    GYH: copied from parm99.dat X -CT-S -X : JCC,7,(1986),230
CT-CT-ss-c3   3    1.00          0.0             3.	    GYH: copied from parm99.dat X -CT-S -X : JCC,7,(1986),230
CT-CT-ss-ss   3    1.00          0.0             3.         GYH: copied from parm99.dat X-CT-S-X : JCC,7,(1986),230
H1-CT-ss-ss   3    1.00          0.0             3.         GYH: copied from parm99.dat X-CT-S-X : JCC,7,(1986)


IMPROPER
CA-CA-CA-HA         1.1          180.0         2.0          Using default value
CA-CA-CA-CA         1.1          180.0         2.0          Using default value
CA-CA-CA-CT         1.1          180.0         2.0          Using default value
CT-CT-C -O          1.1          180.0         2.0          Using default value

NONBON
  CA          1.9080  0.0860             same as ca 
  C           1.9080  0.0860             same as c  
  CT          1.9080  0.1094             same as c3 
  O           1.6612  0.2100             same as o  
  HA          1.4870  0.0157             same as hc 
  HC          1.4870  0.0157             same as hc j
"""
        fout = open('dlm_edit.frcmod','w')
        fout.write(dlm_text)
        fout.close()


	leapscript = ''
        

    # add loadAmberParams commands to the leapscript
	leapscript += 'mods%s = loadAmberParams %s\n'%(p.ResidueName, out_mol2_frcmod) 
	leapscript += 'modsDLM = loadAmberParams dlm_edit.frcmod\n'

        #### NOTE: the choice AMBER forcefield here is hard-coded, alter the code here to change it ####
        ####       Available forcefields will depend on your installation of amber/ambertools !
        #  $ ls -1 $AMBERHOME/dat/leap/cmd/leaprc.ff*
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff02pol.r1
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff02polEP.r1
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff03.r1
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff03ua
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff10
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff99SB     <--
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff99SBildn <---
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff99SBnmr
        #  $AMBERHOME/dat/leap/cmd/leaprc.ff99bsc0
        #  $AMBERHOME/dat/leap/cmd/leaprc.ffAM1
        #  $AMBERHOME/dat/leap/cmd/leaprc.ffPM3

	leapscript += 'source oldff/leaprc.ff99SB\n' 
	leapscript += '\n'

    # add lines to build each residue
	leapscript += tleap_residue_lines(out_mol2, p.ResidueName, p.CNTerm, CappingAtoms, HeadAtom=p.HeadAtom, TailAtom=p.TailAtom)
	leapscript += "# Make a list of info needed for the library\n"
	leapscript += "# This is necessary in order to make sure the parameter modifications are stored to the off file.\n"
	leapscript += "myresidue = { "
	leapscript += "%s mods%s "%(p.ResidueName, p.ResidueName)
	leapscript += "modsDLM "  # Missing bond, angle and torsion parms supplied by DLM and VAV
	leapscript += "}\n"
	leapscript += "saveoff myresidue %s.off\n\n"%p.ResidueName
	leapscript += "quit\n"

        print '### THE leapscript ###'
        print leapscript
        print '### -------------- ###'

        # Write script
        (a, tfile) = tempfile.mkstemp()
        file = open(tfile, 'w')
        file.write(leapscript)
        file.close()

        # Run script
        print commands.getoutput( 'tleap -f %s' % tfile )
        os.remove(tfile)



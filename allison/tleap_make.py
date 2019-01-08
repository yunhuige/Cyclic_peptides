# IMPORTS:{{{
import os, sys, string
import os.path
# }}}

# FUNCTIONS:{{{
def get_sequences(filename):
    """Reads a list of numbered peptide sequences from a text file, formatted like:
    
    # code     sequence
    p0       EEWAREIGAQLRRIADDLNAQYE
    p1       EQWAREIGAQLRXBADXLNAQYE
    ...
    
    and returns a dictionary with entries {p0: 'EEWAREIGAQLRRIADDLNAQYE', ...}
    
    INPUT
    filename - the name of the input file
    
    OUTPUT
    sequences - a dictionary object containing the sequence information.
    """

    # read in the lines from the input file
    fin = open(filename, 'r')
    lines = fin.readlines()
    fin.close()

    sequences = {}
    for line in lines:
        if line[0] != '#':
            fields = line.split()
            sequences[fields[0]] = fields[1]

    return sequences



def make_custom_command_text(command_file, key):
    """Writes custom command portion of tleap script given the command file, the peptide key number, and a peptide sequence string. Saves as text file with format 'pepname_command.txt'"""
    
    fin = open(command_file, 'r')
    command_template = fin.readlines()
    command_lines = [command.replace('pep', key) for command in command_template]
    fin.close()

    fout = open('%s_command.txt'%key, 'w')
    for line in command_lines:
        fout.write(line)
    fout.close()



def write_a_tleap_script(outfile, key, sequence, commands):
    """Writes a tleap script given an output filename, the peptide key number, and a peptide sequence string."""
    
    print 'Writing peptide %s to %s ...'%(key,outfile),

    # Translate the single-letter code in the sequence string to a list of three-letter codes
    three_letter_list = [olc2tlc[sequence[i]] for i in range(len(sequence))]
    
    # Join the three-letter code into a single space-padded string that AMBER can understand
    seq_string = string.joinfields(three_letter_list, ' ')
    
    # Write individualized tleap script to outfile with input file containing necessary loadoff commands called named 'load_customres.txt'
    source_text = """
source leaprc.protein.ff14SB
source leaprc.gaff
gaff = loadamberparams gaff.dat
\n"""

    r = open('load_customres.txt', 'r')
    customres_text = r.read()
    r.close()

    seq_text = """
%s = sequence { %s }
\n"""%(key,seq_string)

    c = open(commands, 'r')
    command_text = c.read()
    c.close()

    wq_text = """
saveAmberParm %s %s.prmtop %s.crd 

quit
\n"""%(key,key,key)

    tleap_text = source_text + customres_text + seq_text + command_text + wq_text
    print 'tleap_text', tleap_text
    
    fout = open(outfile, 'w')
    fout.write(tleap_text)
    fout.close() 

    print '... Done %s.'%key
# }}}

# TRANSLATION TABLE:{{{
# define a translation table from one-letter-code to three-letter code
olc2tlc = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE',
           'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU',
           'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG',
           'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR', 
           'X':'LIG', 'B':'NLE', 'J':'BIP', 'O':'NIB', 'U':'NIP', 
	   'Z':'WCL'}
# }}}

# MAIN PROGRAM:{{{
# Execute program using syntax 'python tleap_make.py filename' where filename=input_file containing peptide codes and sequences

# First, we read in a list of peptide sequences from file
filename = '%s'%(sys.argv[1])
sequences = get_sequences(filename)
print 'sequences', sequences   # sequences = {p0: 'EEWAREIGAQLRRIADDLNAQYE', ...}

# For each sequence, we will write a separate tleap input files to prepare its topology with input file containing custom commands named 'command.txt'. 'command.txt' is deleted to reduce clutter.
for key, sequence in sequences.iteritems():
    outfile = 'tleap_%s.in'%key
    make_custom_command_text('command.txt', key)
    commands = '%s_command.txt'%key
    write_a_tleap_script(outfile, key, sequence, commands)
    os.remove('%s_command.txt'%key)
# }}}

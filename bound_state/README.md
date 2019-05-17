# Wed May 15 10:01:09 EDT 2019

##### Questions:
• Are we going to mutate via rotamers?
• Are we going to use a reference structure to remap and change the chi, phi, psi angles associated with the 3 main res?

### PLAN:
• Use the pdb file: 2axi.pdb 
• Use the sequence list from Table 2. in the paper DOI: 10.1002/cbic.200500452  (only works for natural amino acids)
• Use the rotamers method to mutate the residues according to the sequence provided
• Seq_Table_2.txt is a file that contains only the natural AA sequences. 

# Wed May 15 16:55:05 EDT 2019

### Protocol:

#### Here is the command that you should run (inside this directory):

chimera --nogui --script "mutate_via_rotamers.py 2axi.pdb Seq_Table_2.txt"

### Prospective plans/notes:
• You may also have to change the residue numbers?
• Make sure you have the correct ordering of the paper and sequences in the txt file
• I need to talk to Vince to make sure that I can proceed with the protocol that Yunhui and I have discussed.


# Thu May 16 15:28:14 EDT 2019

### Notes:
• Only thing left for the natural sequences is to manually fix sequence #41 due to the tryptophan being in the position of Phe. (This is the only sequence with this issue.)


# Fri May 17 10:13:53 EDT 2019

### Notes:
Today, we are building a box...









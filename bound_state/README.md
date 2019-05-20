# README.md ——— Lab Notebook 

## Contents:

1. [List of Sequences](#natural-aa-sequences-from-table-2)
2. [Sequence Mapping](#Mapping)


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
`[x]` Talk to Vince to make sure that I can proceed with the protocol that Yunhui and I have discussed.


# Thu May 16 15:28:14 EDT 2019

### Notes:
• Only thing left for the natural sequences is to manually fix sequence #41 due to the tryptophan being in the position of Phe. (This is the only sequence with this issue.)


# Fri May 17 10:13:53 EDT 2019

### Notes:
Today, we are building a box...

Here is the protocol:

```bash

gmx pdb2gmx -f *pdb  -ignh  
gmx editconf -f conf.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p new.top -o solv.gro
gmx grompp -f mdp/minimize.mdp -c solv.gro -p new.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p new.top -pname NA -nname CL -neutral -conc 0.1
gmx make_ndx -f solv_ions.gro

```

# Fri May 17 15:38:20 EDT 2019

• Wrote up a script to reorder the atoms in the structure file (gro) to match the topology (see reorder.py)


# Mon May 20 10:19:01 EDT 2019

#TODO:
• Apply the script to obtain all of the gro files and placed them inside of the reordered the directories according to [the mapping below](#Mapping).

#### Protocol

1. split pdb --> chimera addh to ligand & change DPR to PRO --> ligand.pdb
2. editconf -f ligand.pdb  -o ligand.gro
3. Use reorder.py to get new.gro

#### Commands:

chimera --nogui --script "mutate_via_rotamers.py 2axi.pdb Seq_Table_2.txt"

for i in {0..62};do cd $i;gmx editconf -f ligand.pdb -o ligand.gro;python ../reorder.py ligand.gro /Users/tuc41004/github/Cyclic_peptides/sequences/cb2rr/RUNS/$i/ligand_GMX.top;cd ../;done


For Seq # 41: Change bond angles to match the following:

CD1: -122.32; CG: 67.557

bond angle: CA-CB-CG: 108.0

## Seq # 57 HIS has an error with the number of hydrogens. 

### NOTE: This RUNS directory is incorrect.  Corresponds to the newer mapping.txt file
/Users/tuc41004/github/Cyclic_peptides/protocol_RR/RUNS/0/ligand_GMX.top






------------------------------------


#### Natural AA sequences from Table 2.

| \# | Sequence |
| :--: | :--: |
| 1  | FLWLNKET |
| 6  | FLWLNKEI |
| 7  | FIWLNKET |
| 8  | FIWLNKEI |
| 9  | FIWLNKPT |
| 10 | FIWLNKPI |
| 11 | FIWLKKPI |
| 12 | FLWLGKET |
| 13 | FLWLGKEI |
| 14 | FIWLGKET |
| 15 | FIWLGKEI |
| 16 | FIWLGKPT |
| 17 | FIWLGKPI |
| 18 | FLWLNGET |
| 19 | FLWLKGEI |
| 20 | FIWLKGET |
| 21 | FIWLKGEI |
| 22 | FIWLKGPT |
| 23 | FIWLKGPI |
| 24 | FCWLNKCT |
| 25 | FLWLNYET |
| 26 | FKWLNYET |
| 27 | FEWLNYKT |
| 28 | FKWLNYEF |
| 29 | FKWLNYEY |
| 31 | FKWLNFET |
| 33 | FLWGNYET |
| 34 | FLWNGYET |
| 35 | FLWNPYET |
| 36 | FLWPNYET |
| 37 | FLWVNYET |
| 38 | FKWNGYEF |
| 41 | WKFNGYET |
| 42 | FLWLGKPI |
| 43 | FLWLNYPI |
| 44 | FLWLNFPI |
| 45 | FLWLGFPI |
| 46 | FLWLNKPI |
| 47 | FEWLGKPI |
| 48 | FEWLNWPI |
| 49 | FEWLNFPI |
| 50 | FEWLGFPI |
| 51 | FEWLNYKF |
| 52 | FEWLNFKY |
| 53 | FEWLNFKW |
| 54 | FEWLGFKY |
| 55 | FEWPNFKY |
| 56 | FEWLNWKY |
| 57 | FEWLNHKY |
| 58 | FEWDQFKY |
| 59 | FEWDQFPI |
| 60 | FDWLNFKY |
| 61 | FQWLNFKY |
| 62 | FEWLNFRY |
| 63 | FEWLNFEY |
| 64 | FEWLNFQY |
| 65 | FQWLNFQY |
| 66 | FDWLNFQY |
| 67 | FEWLNWEF |
| 68 | FDWLNFEY |
| 69 | FEWLNWEY |
| 71 | FEWLDWEF |
| 72 | FEWDQWEF |


------------------------------------


### Mapping

Directory # --> Sequence (Table 2.)

```

0 --> 1
1 --> 6
2 --> 7
3 --> 8
4 --> 9
5 --> 10
6 --> 11
7 --> 12
8 --> 13
9 --> 14
10 --> 15
11 --> 16
12 --> 17
13 --> 18
14 --> 19
15 --> 20
16 --> 21
17 --> 22
18 --> 23
19 --> 24
20 --> 25
21 --> 26
22 --> 27
23 --> 28
24 --> 29
25 --> 31
26 --> 33
27 --> 34
28 --> 35
29 --> 36
30 --> 37
31 --> 38
32 --> 41
33 --> 42
34 --> 43
35 --> 44
36 --> 45
37 --> 46
38 --> 47
39 --> 48
40 --> 49
41 --> 50
42 --> 51
43 --> 52
44 --> 53
45 --> 54
46 --> 55
47 --> 56
48 --> 57
49 --> 58
50 --> 59
51 --> 60
52 --> 61
53 --> 62
54 --> 63
55 --> 64
56 --> 65
57 --> 66
58 --> 67
59 --> 68
60 --> 69
61 --> 71
62 --> 72

```



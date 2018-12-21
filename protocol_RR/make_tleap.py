import numpy as np
import sys, os

# Hardcoded Sequences from Table 2 in DOI: 10.1002/cbic.200500452
seq = ['FLWLNKETPP','FLWLNKEIPP','FIWLNKETPP','FIWLNKEIPP',
       'FIWLNKPTPP','FIWLNKPIPP','FIWLKKPIPP','FLWLGKETPP',
       'FLWLGKEIPP','FIWLGKETPP','FIWLGKEIPP','FIWLGKPTPP',
       'FIWLGKPIPP','FLWLNGETPP','FLWLKGEIPP','FIWLKGETPP',
       'FIWLKGEIPP','FIWLKGPTPP','FIWLKGPIPP','FCWLNKCTPP',
       'FLWLNYETPP','FKWLNYETPP','FEWLNYKTPP','FKWLNYEFPP',
       'FKWLNYEYPP','FKWLNFETPP','FLWGNYETPP','FLWNGYETPP',
       'FLWNPYETPP','FLWPNYETPP','FLWVNYETPP','FKWNGYEFPP',
       'WKFNGYETPP','FLWLGKPIPP','FLWLNYPIPP','FLWLNFPIPP',
       'FLWLGFPIPP','FLWLNKPIPP','FEWLGKPIPP','FEWLNWPIPP',
       'FEWLNFPIPP','FEWLGFPIPP','FEWLNYKFPP','FEWLNFKYPP',
       'FEWLNFKWPP','FEWLGFKYPP','FEWPNFKYPP','FEWLNWKYPP',
       'FEWLNHKYPP','FEWDQFKYPP','FEWDQFPIPP','FDWLNFKYPP',
       'FQWLNFKYPP','FEWLNFRYPP','FEWLNFEYPP','FEWLNFQYPP',
       'FQWLNFQYPP','FDWLNFQYPP','FEWLNWEFPP','FDWLNFEYPP',
       'FEWLNWEYPP','FEWLDWEFPP','FEWDQWEFPP']

mimetic = [1,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
        25,26,27,28,29,31,33,34,35,36,37,38,41,42,43,44,45,46,47,
        48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
        67,68,69,71,72]

AA = {'ARG' : 'R','HIS' : 'H','LYS' : 'K','ASP' : 'D',
        'GLU' : 'E','SER' : 'S','THR' : 'T','ASN' : 'N',
        'GLN' : 'Q','CYS' : 'C','SEC' : 'U','GLY' : 'G',
        'PRO' : 'P','ALA' : 'A','VAL' : 'V','ILE' : 'I',
        'LEU' : 'L','MET' : 'M','PHE' : 'F','TYR' : 'Y',
        'TRP' : 'W'}

_seq = []
for i in range(0,len(seq)):
    print seq[i]
    _res = []
    n = 0
    for res in seq[i]:
        for key, val in AA.items():
            if res == val:
                _res.append(key)
                print(_res[n])
                n += 1
    sequence = np.array(_res)
    _seq.append(sequence)

seq = _seq
print(seq)


#fin = open('prototype.txt','r')
fin = open('/Volumes/RMR_4TB/Research/Cyclic_peptide/Cyclic_peptides/protocol/prototype.txt','r')

tleap_text = fin.read()
fin.close()

for i in range(len(seq)):
    s = '%s'%seq[i]
    sequence = s.replace("'","").replace("[","").replace("]","").replace(","," ")
    print(sequence)
    new_tleap = tleap_text.replace('$sequence$',sequence)
    fout = open('tleap_%d.in'%mimetic[i],'w')
    fout.write(new_tleap)
    fout.close()





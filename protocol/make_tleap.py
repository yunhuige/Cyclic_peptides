import numpy as np
import sys, os

seq=['','',...]  # put your sequence here, e.g 'ALA ALA ALA ALA ALA...'

fin = open('prototype.txt','r')
tleap_text = fin.read()
fin.close()

for i in range(len(seq)):
    new_tleap = tleap_text.replace('$sequence$','%s'%seq[i])
    fout = open('tleap_%d.in'%i,'w')
    fout.write(new_tleap)
    fout.close()


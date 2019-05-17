#!/usr/bin/env python

usage="""

python reorder.py [gro] [top]

Description: This script reorders each atom in a structure file (gro) to match the topology.
"""

import sys

if (sys.argv[1:] == ("-h" or "--h" or "--help" or "help")) or (len(sys.argv) != 3):
    print(usage)
    exit(1)

gro = sys.argv[1]
top = sys.argv[2]

newFile = []
with open(top, "r") as file:
    topology = file.readlines()
    # Get the title and the number of atoms from the gro file
    for j in range(0,2):
        structureFile = open(gro, "r")
        sFile = structureFile.readlines()
        newFile.append(str(sFile[j]))
    start, end = 0,0

    # Find the start and end positions inside the topology
    for i in range(len(topology)):
        if topology[i].strip() == "[ atoms ]":
            start = i+2
        if topology[i].strip() == "[ bonds ]":
            end = i-1
    #print(start,end)
    #exit(1)
    # Loop through the gro file to reorder according to the top
    for i in range(start, end):
        ith_line_in_top = list(topology[i].split())
        ndx, res, atom = ith_line_in_top[2], ith_line_in_top[3], ith_line_in_top[4]
        # NOTE: Hardcoded to match the residue number of the top
        RES = "%s"%(str(int(ndx)+20)+res)
        structureFile = open(gro, "r")
        sFile = structureFile.readlines()
        for j in range(2,len(sFile)):
            jth_line_in_struct = list(sFile[j].split())
            if (jth_line_in_struct[0] == RES) and (jth_line_in_struct[1] == atom):
                newFile.append(
                        " "+str(sFile[j]).replace(RES,
                            str(int(jth_line_in_struct[0][:2])-20)+jth_line_in_struct[0][2:]))
                # Extra added space in the line above to match to gromacs spacing style
        # NOTE: Hardcoded for proline to be at the end of the new gro file instead of the beginning
        if ndx == str(10):
            for j in range(2,16):
                jth_line_in_struct = list(sFile[j].split())
                RES = "%s"%(str(int(ndx)+10)+res)
                if (jth_line_in_struct[0] == RES) and (jth_line_in_struct[1] == atom):
                    newFile.append(
                            str(sFile[j]).replace(RES,
                                str(int(jth_line_in_struct[0][:2])-10)+jth_line_in_struct[0][2:]))
    # Get the box size
    newFile.append(str(sFile[-1]))

# Write all the reordered lines from gro to a new file
with open(gro.replace(".gro", "_.gro"), "w") as file:
    for i in range(len(newFile)):
        file.writelines(newFile[i])
    file.close()



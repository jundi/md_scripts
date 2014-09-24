#!/usr/bin/python
#
# An exercise script written in a python course in summer 2014.
# A stupid way to modify .gro file with python. 
#

import sys # sys.argv

file_in=sys.argv[1]
file_out='modified.gro'
coord_trans=[1,2,3]
box_scale=[2,2,2]

# read file to list
lines = [line.rstrip() for line in open(file_in)]


# new coordinates
atoms = []
for line in lines[2:-1]:

    rnum=line[0:5]
    rname=line[5:10]
    aname=line[10:15] 
    anum=line[15:20]
    coord=[line[20:28], line[28:36], line[36:44]]
    if len(line) > 45:
        vel=[line[44:52], line[52:60], line[60:68]]
    else:
        vel=''
    
    # skip waters
    if rname == 'SOL  ':
        continue   
    
    # edit coordinates
    for i in range(3):
        coord[i] = "{0:8.3f}".format(float(coord[i]) + coord_trans[i])
    
    # save modified line
    newline = rnum + rname + aname + anum + ''.join(coord) + ''.join(vel)
    atoms.append(newline)
    


# write new file
with open(file_out, "w+") as new:

    # header and number of atoms
    new.write(lines[0] + '\n')
    new.write(str(len(atoms)) + '\n')
    
    # atoms
    for line in atoms:
        new.write(line + '\n')

    # box coordinates
    box = lines[-1].split()
    newbox = [str(float(box[0])*box_scale[0]), str(float(box[1])*box_scale[1]), str(float(box[2])*box_scale[2])]
    new.write(' '.join(newbox) + '\n')

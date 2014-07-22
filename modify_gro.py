#!/usr/bin/python

import sys # argv

file_in=sys.argv[1]
file_out='modified.gro'
x_trans=1
y_trans=2
z_trans=3
box_x_scale=2
box_y_scale=2
box_z_scale=2



# read file to list
lines = []
with open(file_in, "r") as old:
    for line in old:
        lines.append(line)



# new coordinates
atoms = []
for line in lines[2:-1]:

    rnum=line[0:5]
    rname=line[5:10]
    aname=line[10:15] 
    anum=line[15:20]
    x=line[20:28]
    y=line[28:36]
    z=line[36:44]
    if len(line) > 45:
        vx=line[44:52]
        vy=line[52:60]
        vz=line[60:68]
    

    # skip waters
    if rname == 'SOL  ':
        continue   
    
    
    # edit coordinates
    x_f=float(x)+x_trans
    y_f=float(y)+y_trans
    z_f=float(z)+z_trans
    
    # new coordinates
    x_m="{0:8.3f}".format(x_f)
    y_m="{0:8.3f}".format(y_f)
    z_m="{0:8.3f}".format(z_f)
    
    # save modified line
    if len(line) > 45:
        newline = rnum + rname + aname + anum + x_m + y_m + z_m + vx + vy + vz + '\n'
    else:
        newline = rnum + rname + aname + anum + x_m + y_m + z_m + '\n'

    atoms.append(newline)
    


# write new file
with open(file_out, "w+") as new:


    # header
    for line in lines[:1]:
        new.write(line)


    # number of atoms
    new.write(str(len(atoms)) + '\n')

    
    # atoms
    for line in atoms:
        new.write(line)


    # box coordinates
    line = lines[-1]
    box=line.split()
    box_x = float(box[0])*box_x_scale
    box_y = float(box[1])*box_y_scale
    box_z = float(box[2])*box_z_scale
    newbox = str(box_x) + ' ' + str(box_y) + ' ' + str(box_z) + '\n'
    new.write(newbox)


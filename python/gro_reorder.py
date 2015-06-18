#!/usr/bin/env python3
'''
Reorder atoms in Gromacs coorfinate file (.gro).
'''

import argparse



# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", help='Input (.gro)')
parser.add_argument("-o", help='Output (.gro)')
parser.add_argument("-l", help='File which contains atom numbers in new order (one atom per line)')
parser.add_argument("-r", help='Resname')
args = parser.parse_args()

if args.i:
    input_gro = args.i
if args.o:
    output_gro = args.o
if args.l:
    order_file = args.l
if args.r:
    resname = args.r



# the atoms from old gro-file will written to new file in this order:
order=[]
with open(order_file, 'r') as order_stream:
    for line in order_stream:
        order.append(int(line.split()[0]))



# open coordinate files
oldgro = open(input_gro, 'r')
newgro = open(output_gro, 'w+')



# start reading old gro-file
resnum = -1
atoms = []
line = oldgro.readline()
while line != '':


    # Next residue:
    nextresnum = line[0:5].replace(" ","")
    nextresname = line[5:10].replace(" ","")

    # If this atom starts a new residue, the last residue will be written to
    # file, and the atom list will be cleared.
    if nextresnum != resnum:
        # Check that the amount of atoms is correct
        if len(atoms) == len(order):
            for n in order:
                newgro.write(atoms[n-1])

        # Clear list of atoms
        atoms = []
        # Next residue
        resnum=nextresnum


    # If this atom is part of aresidue which is to be reordered, we add it to
    # list.  If its, not it will written to the file.
    if nextresname == resname and nextresnum.isdigit():
        atoms.append(line)
    else:
        newgro.write(line)

    # Read next line
    line = oldgro.readline()



oldgro.close()
newgro.close()

#!/usr/bin/env python2

import argparse 
import MDAnalysis
import numpy 
import sys


### defaults
resname="LBPA"
#atomnames=["O21'", "C1'", "C3'", "HS'"]
#atomnames=["O21", "C1", "C3", "HS"]
atomnames=["O21'", "C1'", "C3'", "HS'", "O21", "C1", "C3", "HS"]


### Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-c", help='coordinate file', required=True)
parser.add_argument("-a", help='atom list in priority order', nargs='+')
parser.add_argument("-r", help='resname')
args = parser.parse_args()

if args.c:
    coordinatefile = args.c
if args.a:
    if len(args.a) % 4 > 0:
        sys.exit("wrong number of atoms")
    else:
        atomnames = args.a
if args.r:
    resname = args.r


# open file
universe = MDAnalysis.Universe(coordinatefile) 

#number of chiral centers
chirals = len(atomnames)/4

# loop through all residues
resids = universe.selectAtoms("resname " + resname).resids()
for r in resids:

    handedness = []

    for c_itr in range(0,chirals):

        coord = []
        for n in atomnames[c_itr:c_itr+4]:
            atom = universe.selectAtoms("resid " + str(r) + " and name " + n)
            coord.append(atom.positions[0])

        vec12 = coord[1] - coord[0]                 # vector from atom1 to atom2
        vec13 = coord[2] - coord[0]                 # vector from atom1 to atom3
        crs123 = numpy.cross(vec12, vec13)          # cross product
        com123 = (coord[0] + coord[1] + coord[2])/3 # center of mass
        vec4 = com123 - coord[3]                    # vector from com123 to atom4
        dotp = numpy.dot(crs123, vec4)              # dot product

        if dotp < 0:
            handedness.append('S')
        else:
            handedness.append('R')

    print(resname + str(r) + ": " + " ".join(handedness))

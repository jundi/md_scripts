#!/usr/bin/env python2

import argparse 
import MDAnalysis
import numpy 
import sys


### DEFAULTS
coordinatefile="confout.gro"
resname="LBPA"
atomnames=["O21'", "C1'", "C3'", "HS'", "O21", "C1", "C3", "HS"]



### PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("-c", help='Coordinate file')
parser.add_argument("-a", help='Atoms bound to the chiral center (in order of'
                               ' priority). If there is n chiral centers, 4n '
                               'atoms should be given.', nargs='+')
parser.add_argument("-r", help='Resname of the molecule to be checked.')
args = parser.parse_args()

if args.c:
    coordinatefile = args.c
if args.a:
    if len(args.a) % 4 > 0:
        sys.exit("Wrong number of atoms.")
    else:
        atomnames = args.a
if args.r:
    resname = args.r



### OPEN FILE AND CHECK CHIRALS
universe = MDAnalysis.Universe(coordinatefile) 

# number of chiral centers
chirals = len(atomnames)/4

# loop through all residues
resids = universe.selectAtoms("resname " + resname).resids()
for r in resids:

    handedness = [] # handedness of each chiral centers in residue

    # loop through all chiral centers
    for c_itr in range(0,chirals):

        # get coordinates of four atoms
        coord = []
        for n in atomnames[c_itr:c_itr+4]:
            atom = universe.selectAtoms("resid " + str(r) + " and name " + n)
            coord.append(atom.positions[0])

        # get the normal of plane spanned by three atoms
        vec12 = coord[1] - coord[0]                 # vector from atom1 to atom2
        vec13 = coord[2] - coord[0]                 # vector from atom1 to atom3
        crs123 = numpy.cross(vec12, vec13)          # plane normal vector

        # which side of the plane is the fourth atom?
        com123 = (coord[0] + coord[1] + coord[2])/3 # center of mass
        vec4 = com123 - coord[3]                    # vector from the plane to atom4

        # Is the the plane normal vector parallel or antiparallel with vec4?
        dotp = numpy.dot(crs123, vec4)              # dot product
        if dotp < 0:
            handedness.append('S')
        else:
            handedness.append('R')

    print(resname + str(r) + ": " + " ".join(handedness))

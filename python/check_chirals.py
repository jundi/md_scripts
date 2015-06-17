#!/usr/bin/env python2
"""
Checks the R/S-handedness of molecules stereocenters based on the coordinates of
the four atoms bound to a chiral center.
"""
import argparse 
import MDAnalysis
import numpy 
import sys



### DEFAULTS
coordinatefile="confout.gro"
resname="LBPA"
atomnames=["O21", "C1", "C3", "HS", "O21'", "C1'", "C3'", "HS'"]



### PARSE ARGUMENTS
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", help='Coordinate file')
parser.add_argument("-a", help='Atoms bound to the chiral center (in order of'
                               ' priority). If there is n chiral centers, 4n '
                               'atoms should be given.', nargs='+')
parser.add_argument("-r", help='Resname of the molecule to be checked.')
args = parser.parse_args()

if args.c:
    coordinatefile = args.c
else:
    print("Using default coordinate file: " + coordinatefile)
if args.a:
    if len(args.a) % 4 > 0:
        sys.exit("Wrong number of atoms.")
    else:
        atomnames = args.a
else:
    print("Using default atoms: " + resname)
if args.r:
    resname = args.r
else:
    print("Using default resname: " + ", ".join(atomnames))



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
        for n in atomnames[4*c_itr:4*c_itr+4]:
            atom = universe.selectAtoms("resid " + str(r) + " and name " + n)
            coord.append(atom.positions[0])

        # get the normal of plane spanned by three atoms
        vec12 = coord[1] - coord[0]                 # vector from atom1 to atom2
        vec13 = coord[2] - coord[0]                 # vector from atom1 to atom3
        crs123 = numpy.cross(vec12, vec13)          # plane normal vector

        # which side of the plane is the fourth atom?
        com123 = (coord[0] + coord[1] + coord[2])/3 # center of mass
        vec4 = coord[3] - com123                    # vector from the plane to atom4

        # Is the the plane normal vector parallel or antiparallel with vec4?
        dotp = numpy.dot(crs123, vec4)              # dot product
        if dotp < 0:
            # angle bigger than 90 degrees -> antiparallel -> "Sinister"
            handedness.append('S')
        else:
            # angle smaller than 90 degrees -> parallel -> "Rectus"
            handedness.append('R')

    print("{0}{1:5}{2}".format(resname, str(r), " ".join(handedness)))

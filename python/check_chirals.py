#!/usr/bin/env python2
"""
Checks the R/S-handedness of molecules stereocenters based on the coordinates of
the four atoms bound to a chiral center.
"""
import argparse 
import MDAnalysis
import numpy 
import sys



### ATOMDICT
# some atom lists hard-coded
atomdict={}
atomdict["LBPA22"] = ["O21", "C1", "C3", "HS", "O21'", "C1'", "C3'", "HS'"]
atomdict["LBPA33"] = ["OC2", "C1", "C3", "H2A", "OC2'", "C1'", "C3'", "H2A'"]



### DEFAULTS
coordinatefile="confout.gro"
resname="LBPA"
atomnames=atomdict["LBPA22"]



### PARSE ARGUMENTS
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", help='Coordinate file')
parser.add_argument("-a", help='''Atoms bound to the chiral center
                                  (in order of priority). If there is n chiral
                                  centers, 4n atoms should be given.  If only
                                  one argument is given, some of the hard coded
                                  lists are used.''', nargs='+')
parser.add_argument("-r", help='Resname of the molecule to be checked.')
parser.add_argument("-l", help='Show hard coded atom lists.', action='store_true')
args = parser.parse_args()

if args.l:
    for i in atomdict:
        print("{0}: {1}".format(i, " ".join(atomdict[i])))
    sys.exit()

if args.c:
    coordinatefile = args.c
else:
    print("Using default coordinate file: " + coordinatefile)

if args.r:
    resname = args.r
else:
    print("Using default resname: " + resname)

if args.a:
    if len(args.a) == 1:
        if args.a[0] in atomdict:
            atomnames = atomdict[args.a[0]]
            print("Using atoms: " + ", ".join(atomnames))
        else:
            sys.exit("Atom list {0} is not known.".format(args.a[0]))

    elif len(args.a) % 4 > 0:
        sys.exit("Wrong number of atoms.")
    else:
        atomnames = args.a
else:
    print("Using default atoms: " + ", ".join(atomnames))



### OPEN FILE AND CHECK CHIRALS
universe = MDAnalysis.Universe(coordinatefile) 

# number of chiral centers
chirals = len(atomnames)/4

# loop through all residues
residues = universe.selectAtoms("resname " + resname).residues
for r in residues:

    handedness = [] # handedness of each chiral centers in residue

    # loop through all chiral centers
    for c in range(0,chirals):

        # get coordinates of four atoms
        coord = []
        for n in atomnames[4*c:4*c+4]:
            atom = getattr(r, n)
            coord.append(atom.position)

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

    print("{0}{1:5}{2}".format(resname, str(r.id), " ".join(handedness)))
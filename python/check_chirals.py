#!/usr/bin/env python2
"""
Checks the R/S-handedness of molecules stereocenters based on the coordinates
of the chiral center and four atoms bound to it. This script can not handle
PBC, so take care that molecules in your coordinate files are not broken.
"""
import argparse 
import MDAnalysis
import numpy 
import sys
import os.path



### ATOMDICT
# some default atom lists hard coded
atomdict={}
# LBPA 2-2' isoform
atomdict["LBPA22"] = ["C2", "O21", "C1", "C3", "HS", "C2'", "O21'", "C1'", "C3'", "HS'"]
# LBPA 3-3' isoform
atomdict["LBPA33"] = ["C2", "OC2", "C1", "C3", "H2A", "C2'", "OC2'", "C1'", "C3'", "H2A'"]



### PARSE ARGUMENTS
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", help='''Coordinate file. Any format supported by
        MDAnalysis (http://www.mdanalysis.org/) should work.''')
parser.add_argument("-a", help='''Names of chiral center atom and atoms bound
        to the chiral center (in order of priority).  If there is n chiral
        centers, 4n atoms should be given. If only one argument is given, some
        of the hard code lists are used.''', nargs='+')
parser.add_argument("-r", help='Resname of the molecule to be checked.')
parser.add_argument("-l", help='Show hard coded atom lists.', action='store_true')
parser.add_argument("-s", help='''Set configuration to R or S. This option
        takes as many arguments as there is chiral centers in the molecule.''',
nargs='+', choices=['R', 'S'])
parser.add_argument("-o", help='Output file (if using -s)')
args = parser.parse_args()

# list hard coded atom lists
if args.l:
    for i in atomdict:
        print("{0}: {1}".format(i, " ".join(atomdict[i])))
    sys.exit()

# set coordinate file
if args.c:
    coordinatefile = args.c
else:
    coordinatefile="confout.gro"
    print("Using default coordinate file: " + coordinatefile)

# set resname
if args.r:
    resname = args.r
else:
    resname="LBPA"
    print("Using default resname: " + resname)

# set atom list
if args.a:
    if len(args.a) == 1:
        if args.a[0] in atomdict:
            atomnames = atomdict[args.a[0]]
            print("Using atoms: " + ", ".join(atomnames))
        else:
            sys.exit("error: Atom list {0} is not known.".format(args.a[0]))

    elif len(args.a) % 5 > 0:
        sys.exit("error: Wrong number of atoms.")
    else:
        atomnames = args.a
else:
    atomnames=atomdict["LBPA22"]
    print("Using default atoms: " + ", ".join(atomnames))

# set list of labels
if args.s:
    if len(args.s) != (len(atomnames) / 5):
        sys.exit("error: The number of arguments for -s should be equal to number of chiral centers in the molecule.")
    else:
        label = args.s
    if not (args.o):
        sys.exit("error: No output file given.")


# set output file
if args.o:
    if os.path.exists(args.o):
        exit("error: File {} already exists".format(args.o))
    else:
        outputfile = args.o



### OPEN FILE AND CHECK CHIRALS
universe = MDAnalysis.Universe(coordinatefile) 

# number of chiral centers
chirals = len(atomnames)/5

# loop through all residues
residues = universe.selectAtoms("resname " + resname).residues
for r in residues:

    handedness = [] # handedness of each chiral centers in residue

    # loop through all chiral centers
    for c in range(0,chirals):

        # get coordinates of four atoms
        atoms = MDAnalysis.core.AtomGroup.AtomGroup([])
        for n in atomnames[5*c:5*c+5]:
            atoms = atoms + getattr(r, n)
        coord = atoms.positions

        # get the normal of plane spanned by three atoms
        vec12 = coord[2] - coord[1]                 # vector from atom1 to atom2
        vec13 = coord[3] - coord[1]                 # vector from atom1 to atom3
        crs123 = numpy.cross(vec12, vec13)          # plane normal vector

        # which side of the plane is the fourth atom?
        com123 = (coord[1] + coord[2] + coord[3])/3 # center of mass
        vec4 = coord[4] - com123                    # vector from the plane to atom4
        dotp = numpy.dot(crs123, vec4)              # dot product
        if dotp < 0:
            # angle bigger than 90 degrees -> antiparallel -> "Sinister"
            handedness.append('S')
        else:
            # angle smaller than 90 degrees -> parallel -> "Rectus"
            handedness.append('R')

        # Switch handedness?
        if args.s and len(label) > c and (handedness[-1] != label[c]):
            # move fourth atom
            atoms[4].position = coord[4] - 2*vec4           
            # move chiral center
            atoms[0].position = coord[0] - 2*(coord[0] - com123)

            print("Switching {0}. chiral center in {1} to {2}".format(str(c+1), resname+str(r.id), label[c]))

    if not (args.s):
        print("{0}{1:5}{2}".format(resname, str(r.id), " ".join(handedness)))

# write new file with modified chirals
if args.s:
    universe.atoms.write(outputfile)

#!/usr/bin/env python2
"""
Swaps coordinates of two residues.
"""

import MDAnalysis
import argparse 

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-f", help='''Input coordinate file''')
parser.add_argument("-o", help='''Output coordinate file''')
parser.add_argument("-r1", help='''Residue number 1''')
parser.add_argument("-r2", help='''Residue number 2''', nargs='+')
args = parser.parse_args()

infile=args.f
resnum1=args.r1
resnums2=args.r2
if args.o:
    outfile=args.o
else:
    outfile='swap'


for resnum2 in resnums2:
    # open file
    u = MDAnalysis.Universe(infile)
    
    # get positions
    residue1=u.select_atoms("resid {}".format(resnum1))
    residue2=u.select_atoms("resid {}".format(resnum2))
    positions1=residue1.positions
    positions2=residue2.positions
    
    # change positions
    residue1.set_positions(positions2)
    residue2.set_positions(positions1)
    
    # write file
    s=u.select_atoms("segid SYSTEM")
    outfilename="{}_{}_to_{}.gro".format(outfile,str(resnum2), str(resnum1))
    s.write(outfilename)

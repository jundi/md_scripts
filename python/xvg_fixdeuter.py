#!/usr/bin/env python

'''Combines saturated and unsaturated deuter.xvg -plots, and fixes atom
numbering.

For example POPC oleyl tail:
xvg_fixdeuter -f saturated.xvg -u unsaturated.xvg -a 9 10 -o fixed.xvg
'''


import argparse
import xvgio

# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', required=True)
parser.add_argument('-u' )
parser.add_argument('-a', nargs='+')
parser.add_argument('-o', required=True)
args = parser.parse_args()

filename = args.f
unsat_filename = args.u
unsat_atoms = args.a
output_filename = args.o

# read file
tail,comment=xvgio.read(filename)

# fix numbering
tail[:,0]=tail[:,0]+1

if args.u:
    # read file
    unsat,comment2=xvgio.read(unsat_filename)
    # fix numbering
    unsat[:,0]=unsat[:,0]+1
    # merge saturated and unsaturated
    for atom in unsat_atoms:
        index=int(atom)-2
        tail[index,1]=unsat[index,1]

# write file
xvgio.write(output_filename,tail,comment)

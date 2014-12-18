#!/bin/env python2
#
# Reads two .xvg-plots returns difference of their y-axis
# 
# Requires GromacsWrapper and argparse.
#

import argparse
import gromacs.formats


### Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f1", help='file1')
parser.add_argument("-f2", help='file2')
parser.add_argument("-o", help='output')
args = parser.parse_args()

if args.f1:
  filename1 = args.f1
if args.f2:
  filename2 = args.f2
if args.o:
  outfilename = args.o


### Calculate difference
file1 = gromacs.formats.XVG(filename1).array
file2 = gromacs.formats.XVG(filename2).array
diff = file1 
diff[1:,:] = file2[1:,:]-file1[1:,:]


### Write output
outfile = gromacs.formats.XVG()
outfile.set(diff)
outfile.write(outfilename)

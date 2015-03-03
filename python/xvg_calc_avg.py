#!/usr/bin/env python2
#
# Calculate average and standard deviation of xvg-plot.
# 
# Requires GromacsWrapper and argparse.
#

import argparse
import gromacs.formats
import numpy


### Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help='input', required=True)
parser.add_argument("-b", help='first frame (ps)')
args = parser.parse_args()
if args.f:
  filename = args.f
if args.b:
  firstfr = float(args.b)
else:
  firstfr = 0


# get file
xvg = gromacs.formats.XVG(filename).array

# skip frames and remove time column
xvg = xvg[1:,xvg[0] >= firstfr]

# calculate mean and avg for all columns
for column in xvg:
    print(str(column.mean()) + ' ' + str(column.std()))


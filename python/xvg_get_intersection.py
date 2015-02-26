#!/usr/bin/env python2
#
# Calculate intersection of two first columns in xvg-file. If many files are
# given, their average and standard deviation are calculated.
# 
# Requires GromacsWrapper and argparse.
#

import argparse
import gromacs.formats
import numpy


### Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help='input', required=True, nargs='*')
parser.add_argument("-p", help='Predicted value (in case of multiple intersections)')
args = parser.parse_args()
if args.f:
  filename = args.f
if args.p:
  predicted = float(args.p)
else:
  predicted = 1


### Define functions
def get_zerocrossings(filename,predicted):
    # get file
    xvgfile = gromacs.formats.XVG(filename).array

    # interpolate
    x=numpy.arange(xvgfile[0,0],xvgfile[0,-1],0.0001)
    y1=numpy.interp(x,xvgfile[0,:],xvgfile[1,:])
    y2=numpy.interp(x,xvgfile[0,:],xvgfile[2,:])
    
    # find intersection
    y=y1-y2
    zero_crossings = x[numpy.where(numpy.diff(numpy.sign(y)))[0]]
    chosen_index=numpy.argmin(abs(zero_crossings-predicted))
    return zero_crossings[chosen_index]


### Main program
intersects=[]
for f in filename:
    intersects.append(get_zerocrossings(f,predicted))
intersects=numpy.array(intersects)

print(intersects.mean())
print(intersects.std())

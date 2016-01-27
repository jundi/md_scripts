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
parser.add_argument("-f", help='input', required=False, nargs='*')
parser.add_argument("-f1", help='input', required=False, nargs='*')
parser.add_argument("-f2", help='input', required=False, nargs='*')
parser.add_argument("-p", help='Predicted value (in case of multiple intersections)')
args = parser.parse_args()
if args.f:
    filename = args.f
    if args.f1:
        print("error1")
        exit()
    if args.f2:
        print("error1")
        exit()
if args.p:
    predicted = float(args.p)
else:
    predicted = 1
if args.f1:
    filename1 = args.f1
if args.f2:
    filename2 = args.f2

if not (len(filename1) == len(filename2)):
    print("error2")
    exit()



### Define functions
def get_intersections(filename1,predicted,filename2):
    if filename2:
        # get files
        xvgfile1 = gromacs.formats.XVG(filename1).array
        xvgfile2 = gromacs.formats.XVG(filename2).array
        # interpolate
        x=numpy.arange(xvgfile1[0,0],xvgfile1[0,-1],0.0001)
        y1=numpy.interp(x,xvgfile1[0,:],xvgfile1[1,:])
        y2=numpy.interp(x,xvgfile2[0,:],xvgfile2[1,:])
    else:
        # get file
        xvgfile = gromacs.formats.XVG(filename1).array
        # interpolate
        x=numpy.arange(xvgfile[0,0],xvgfile[0,-1],0.0001)
        y1=numpy.interp(x,xvgfile[0,:],xvgfile[1,:])
        y2=numpy.interp(x,xvgfile[0,:],xvgfile[2,:])
    
    # find intersection
    y=y1-y2
    intersections = x[numpy.where(numpy.diff(numpy.sign(y)))[0]]
    intersections = intersections[intersections != x[0]]
    intersections = intersections[intersections != x[-1]]
    # Choose the intersection you were looking for
    intersection = intersections[numpy.argmin(abs(intersections-predicted))]
    print(filename1 + ':    ' + str(intersection) + '   ' + str(intersections))
    return intersection


### Main program
intersects=[]
if args.f:
    for f in filename:
        intersects.append(get_intersections(f,predicted))
else: 
    for f1, f2 in zip(filename1, filename2):
        intersects.append(get_intersections(f1,predicted,f2))
intersects=numpy.array(intersects)

print('')
print('Average: ' + str(intersects.mean()))
print('Std: ' + str(intersects.std()))

#!/bin/python
#
"""
Transfers the x-intercept (or y-intercept) of xvg-plot to a chosen data point.
The default is the last point of dataset (index=-1).
"""
import argparse
import numpy
import xvgio


### arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f',help="Input file")
parser.add_argument('-o',help="Output file")
parser.add_argument('-a',help="Axis",type=int,default=1)
parser.add_argument('-n',help="Index of point to be set to zero",type=int,default=-1)


### parse arguments
args = parser.parse_args()
infile_name = args.f 
if args.o:
    outfile_name = args.o
else:
    outfile_name = '.'.join(infile_name.split('.')[:-1]) + '_trans.' + infile_name.split('.')[-1]
axis = args.a
index = args.n


### read file
data, comment = xvgio.read(infile_name)


### edit data
data[:,axis]=data[:,axis]-data[index,axis]


### write
xvgio.write(outfile_name, data, comment)

#!/bin/python
"""
Cut part of xvg-files. 
"""

import argparse
import numpy
import xvgio


### parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f',help="Input file")
parser.add_argument('-o',help="Output file")
parser.add_argument('-b',help="First data point",type=float)
parser.add_argument('-e',help="Last data point",type=float)

args = parser.parse_args()
infile_name = args.f 
outfile_name = args.o
first = args.b
last = args.e

# read file
data, comment = xvgio.read(infile_name)

# cut data
indices = numpy.all([data[:,0]>=first, data[:,0]<=last], axis=0)
cut_data = data[indices,:]

# write file
xvgio.write(outfile_name, cut_data, comment)

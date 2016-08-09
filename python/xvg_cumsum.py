#!/bin/python
# 
# Cumulative sum of columns.
#

import argparse
import numpy
import xvgio


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='+')
parser.add_argument('-o')

args = parser.parse_args()
infile_names = args.f 
outfile_name = args.o

# read file
data, comment = xvgio.read(infile_name)

data=numpy.cumsum(data,axis=0)

### write
xvgio.write(outfile_name, data, comment)

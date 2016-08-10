#!/bin/python
""" 
Cumulative sum of columns. Columns to be summed can be chosen with -c.
Default is all except first column.
"""

import argparse
import numpy
import xvgio


### parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f',help="Input file")
parser.add_argument('-o',help="Output file")
parser.add_argument('-c',help="Columns",nargs='+',type=int)

args = parser.parse_args()
infile_name = args.f 
outfile_name = args.o
columns = args.c

# read file
data, comment = xvgio.read(infile_name)

if not (columns):
    columns=range(1, data.shape[1])

result=data # init array
for c in columns:
    result[:,c]=numpy.cumsum(data[:,c],axis=0)

### write
xvgio.write(outfile_name, result, comment)

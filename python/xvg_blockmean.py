#!/bin/python
#
# Computes mean of blocks of width N.
#

import argparse
import sys
import os
import xvgio
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True)
parser.add_argument('-o')
parser.add_argument('-n', type=int, default=10)

args = parser.parse_args()

width = int(args.n)
if width < 2:
    sys.exit('n < 2')

infile_name = args.f 

if args.o:
    outfile_name = args.o
else:
    split = os.path.splitext(infile_name)
    outfile_name = split[0] + '_' + str(width) + '-' + 'mean' + split[1]



### read file
data, comment = xvgio.read(infile_name)

### compute means
rows = data.shape[0]
columns = data.shape[1]
newdata = np.zeros([int(rows/width-1), columns])
for r in range(int(rows/width-1)):
    for c in range(columns):
        newdata[r,c] = np.mean(data[r*width:r*width+width-1,c])



### write
xvgio.write(outfile_name, newdata, comment)

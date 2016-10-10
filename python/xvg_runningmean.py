#!/bin/python
#
# Computes simple running mean.
#

import argparse
import sys
import math
import os
import xvgio
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True)
parser.add_argument('-o')
parser.add_argument('-n', type=int, default=9)

args = parser.parse_args()

width = int(args.n)
if width % 2 != 1:
    sys.exit('n has to be odd integer')

infile_name = args.f 

if args.o:
    outfile_name = args.o
else:
    split = os.path.splitext(infile_name)
    outfile_name = split[0] + '_' + str(width) + '-' + 'mean' + split[1]



### read file
data, comment = xvgio.read(infile_name)


### compute means
columns = data.shape[1]
newdata=data[math.floor(width/2):-math.floor(width/2),:]
for c in range(1,columns):
    column = np.cumsum(data[:,c])
    column[width:] = column[width:] - column[:-width]
    column = column[width-1:] / width
    newdata[:,c]=column


### write
xvgio.write(outfile_name, newdata, comment)

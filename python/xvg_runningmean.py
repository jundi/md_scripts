#!/bin/python
#
# Computes running mean
#

import argparse
import sys
import os
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
new_lines = []
data = []
with open(infile_name,'r') as infile:
    for line in infile:
        if line.startswith('#') or line.startswith('@'):
            new_lines.append(line)
        else:
            points = [float(x) for x in line.split()]
            data.append(points)


### compute means
data = np.array(data) 
rows = data.shape[0]
columns = data.shape[1]
newdata = np.zeros([int(rows/width-1), columns])
for r in range(int(rows/width-1)):
    for c in range(columns):
        newdata[r,c] = np.mean(data[r*width:r*width+width-1,c])



### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in newdata:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

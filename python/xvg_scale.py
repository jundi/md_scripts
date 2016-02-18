#!/bin/python
#
# Reads scales xvg-plots
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f')
parser.add_argument('-o')
parser.add_argument('-x', default=1)
parser.add_argument('-y', default=1)

args = parser.parse_args()
infile_name = args.f 
outfile_name = args.o
xscale = args.x
yscale = args.y



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


### scale
data = np.array(data) 
x = float(xscale)
y = float(yscale)
data[:,1:] = y*data[:,1:]
data[:,0] = x*data[:,0]


### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

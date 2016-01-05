#!/bin/python
#
# Move data on axis.
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True)
parser.add_argument('-o')
parser.add_argument('-a', help="axis")
parser.add_argument('-n', help="points")

args = parser.parse_args()
infile_name = args.f 
if args.o:
    outfile_name = args.o
else:
    outfile_name = '.'.join(infile_name.split('.')[:-1]) + '_moved' + '.' + infile_name.split('.')[-1]
if args.n:
    n = int(args.n)
else:
    n = 1000
if args.a:
    axis = int(args.a)
else:
    axis = 0


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


### move data
data = np.array(data) 
data[:,axis] = data[:,axis] + n


### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

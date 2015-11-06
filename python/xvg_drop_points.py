#!/bin/python
#
# Reads scales xvg-plots
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True)
parser.add_argument('-o')
parser.add_argument('-c', help="column")
parser.add_argument('-n', help="every nth point will be saved")

args = parser.parse_args()
infile_name = args.f 
if args.o:
    outfile_name = args.o
else:
    outfile_name = '.'.join(infile_name.split('.')[:-1]) + '_droppoints' + '.' + infile_name.split('.')[-1]
if args.c:
    col = int(args.c)
else:
    col = 2
if args.n:
    n = int(args.n)
else:
    n = 10



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


### set some points to zero
data = np.array(data) 
for i in range(n-1):
    data[i::n,col] = 0


### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

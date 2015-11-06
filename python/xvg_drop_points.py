#!/bin/python
#
# Get datapoints where remainder of x divided by n is zero (x%n == 0).
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True)
parser.add_argument('-o')
parser.add_argument('-n', help="divisor")

args = parser.parse_args()
infile_name = args.f 
if args.o:
    outfile_name = args.o
else:
    outfile_name = '.'.join(infile_name.split('.')[:-1]) + '_droppoints' + '.' + infile_name.split('.')[-1]
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


### find points to be printed
data = np.array(data) 
new_data = data[data[:,0] % n == 0]


### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in new_data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

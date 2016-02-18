#!/bin/python
#
# Get datapoints where remainder of x divided by n is zero (x%n == 0).
#

import argparse
import sys
import numpy as np
import decimal as dec

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
    n = float(args.n)
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
            remainder = abs(points[0]%n)
            if ((remainder < 1E-12) or (abs(remainder-n) < 1E-12)):
                data.append(points)


### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for i in row:
            line=line + str('%11e' % i) + ' '
        outfile.write(line + '\n')


#!/bin/python
#
# Scale y so that integral from b to e is 1.
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f',help='xvg-file',required=True)
parser.add_argument('-o',help='output file',required=True)
parser.add_argument('-b',help='begin',type=float,required=True)
parser.add_argument('-e',help='end',type=float,required=True)
parser.add_argument('-v',help='Be verbose.',action='store_true')

args = parser.parse_args()
infile_name = args.f 
args = parser.parse_args()
outfile_name = args.o 
args = parser.parse_args()
begin = args.b 
args = parser.parse_args()
end = args.e 

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

### integrate
data = np.array(data) 
x = data[(data[:,0] > begin) & (data[:,0] < end), 0]
y = data[(data[:,0] > begin) & (data[:,0] < end), 1]
area = np.trapz(y,x)
if args.v:
    print('scale factor: ' + str(area))

### scale
data[:,1:] = data[:,1:]/area

### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

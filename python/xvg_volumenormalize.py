#!/bin/python
#
# Divides y1 and y2 with 4/3*pi*x^3.
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f')
parser.add_argument('-o')

args = parser.parse_args()
infile_name = args.f 
outfile_name = args.o



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


### compute
data = np.array(data) 
binwidth = data[1,0] - data[0,0]
print("Detected bin = ", + binwidth)
pi = np.pi

# volume from x
vol = 4/3*pi * ( np.power(data[:,0]+binwidth,3) - np.power(data[:,0],3) )

# y1
data[:,1] = data[:,1]/vol[:]
# y2
data[:,2] = data[:,2]/vol[:]



### write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

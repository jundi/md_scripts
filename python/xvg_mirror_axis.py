#!/bin/python
#
# Reads .xvg-plot and mirrors x-axis
#

import sys
import numpy as np

# files
infile_name = sys.argv[1]
outfile_name = '.'.join(infile_name.split('.')[:-1]) + '-mirror.' + infile_name.split('.')[-1]


# read file
new_lines = []
data = []
with open(infile_name,'r') as infile:
    for line in infile:
        if line.startswith('#') or line.startswith('@'):
           new_lines.append(line)
        else:
           points = [float(x) for x in line.split()]
           data.append(points)


# transform
data = np.array(data) 
data[:,0] = -data[:,0]
data=data[::-1]


# write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in data:
        line=''
        for n in row:
            line=line + str('%11e' % n) + ' '

        outfile.write(line + '\n')

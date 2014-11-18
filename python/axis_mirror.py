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
x=[]
y=[]
with open(infile_name,'r') as infile:
    for line in infile:
        if line.startswith('#') or line.startswith('@'):
           new_lines.append(line)
        else:
           data_point=line.split()
           x.append(float(data_point[0]))
           y.append(float(data_point[1]))


# transform
x = np.array(x) 
x = -x[::-1]
y = np.array(y) 
y = -y[::-1]


# print
xy = np.vstack([x,y]).transpose()
print(xy)


# write
with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in xy:
        x='%11e' % row[0]
        y='%11e' % row[1]
        line=str(x) + ' ' + str(y)
        outfile.write(line + '\n')

#!/bin/python
#
# Reads .xvg-plot and sets the y-coordinate of the last data point to zero.
#

import sys
import numpy as np

infile_name = sys.argv[1]
outfile_name = '.'.join(infile_name.split('.')[:-1]) + '-trans.' + infile_name.split('.')[-1]

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


x = np.array(x) 
y = np.array(y) 
y = y-y[-1]
xy = np.vstack([x,y]).transpose()
print(xy)

with open(outfile_name,'w') as outfile:
    for line in new_lines:
        outfile.write(line)
    for row in xy:
        x='%11e' % row[0]
        y='%11e' % row[1]
        line=str(x) + ' ' + str(y)
        outfile.write(line + '\n')

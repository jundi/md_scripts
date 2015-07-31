#!/bin/python
"""Reads .xvg-plot and mirrors x-axis"""

import sys
import argparse
import numpy as np


# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', required=True, nargs='+')
args = parser.parse_args()
infile_names = args.f


for infile_name in infile_names:
    
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
    outfile_name = '.'.join(infile_name.split('.')[:-1]) + '-mirror.' + infile_name.split('.')[-1]
    with open(outfile_name,'w') as outfile:
        for line in new_lines:
            outfile.write(line)
        for row in data:
            line=''
            for n in row:
                line=line + str('%11e' % n) + ' '
    
            outfile.write(line + '\n')

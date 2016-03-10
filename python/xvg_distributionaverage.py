#!/bin/python
#
# Computes mean of distribution.
#

import argparse
import sys
import numpy as np


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='*',help='xvg-file')
parser.add_argument('-v',help='Be verbose.',action='store_true')

args = parser.parse_args()
file_names = args.f 


averages=np.array([])
yvalues=np.array([])
for file_name in file_names:
    ### read file
    new_lines = []
    data = []
    with open(file_name,'r') as infile:
        for line in infile:
            if line.startswith('#') or line.startswith('@'):
               new_lines.append(line)
            else:
               points = [float(x) for x in line.split()]
               data.append(points)
    
    
    ### compute
    data = np.array(data) 
    avg = sum(data[:,0]*data[:,1])/sum(data[:,1])
    if args.v:
        print(avg)
    averages = np.append(averages, avg)

    ### y-value at average
    yvalue = data[np.argmin(abs(data[:,0]-avg)), 1]
    if args.v:
        print(yvalue)
    yvalues = np.append(yvalues, yvalue)



### print
print('mean: ' + str(averages.mean()) + ' +- ' + str(averages.std()))
print('y-value: ' + str(yvalues.mean()) + ' +- ' + str(yvalues.std()))

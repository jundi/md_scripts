#!/bin/python
"""
Concatenate multiple xvg-files. 
"""

import argparse
import numpy
import xvgio


### parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', nargs='+',help="Input files")
parser.add_argument('-o',help="Output file")

args = parser.parse_args()
infile_names = args.f 
outfile_name = args.o

datasets = []
for infile_name in infile_names:
    # read file
    data, comment = xvgio.read(infile_name)
    datasets.append(data)


# concatenate
alldata = numpy.empty((0,numpy.size(datasets[0][0])))
for dataset in datasets:
    alldata=numpy.concatenate((alldata,dataset),axis=0)


# renumber first column
binwidth = alldata[1,0] - alldata[0,0]
print("Detected bin = ", + binwidth)
n=0
for i in alldata:
    i[0]=n
    n=n+binwidth


### write
xvgio.write(outfile_name, alldata, comment)

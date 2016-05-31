#!/bin/python
# 

import argparse
import numpy
import xvgio


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='+')
parser.add_argument('-b', nargs='+', type=int)
parser.add_argument('-o')

args = parser.parse_args()
infile_names = args.f 
first_frames = args.b 
outfile_name = args.o

datasets = []
for infile_name, first_frame in zip(infile_names,first_frames):
    # read file
    data, comment = xvgio.read(infile_name)
    # save part of data
    data = data[data[:,0]>=500000]
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

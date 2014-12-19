#!/usr/bin/python
#
# This script copies plots from multiple .xvg files into one file. The legends
# for plots are concluded from paths of the files. So if your directory has
# contents:
#
#   ./POPC/charge.xvg
#   ./CHOL/charge.xvg
#   ./Protein/charge.xvg
#
# and you merge plots from these files, their legends will be "POPC", "CHOL"
# and "Protein".
#

import argparse
import os
import numpy as np



#--- args ---

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='+')
parser.add_argument('-o')
args = parser.parse_args()


# output file
outfile = args.o
print(outfile)

# list of files
filepaths = args.f
print(filepaths)




#--- FUNCTIONS ---

def read_xvg(filename):
    '''Reads one xy-format .xvg -file.'''
    # Open file
    x=[]
    y=[]
    with open(filename, 'r') as xvg:
        for line in xvg:

            # Skip comment/header lines 
            if (line.startswith('#') or  line.startswith('@')):
                continue

            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))
     
    return  x,y



#--- MAIN ---

# list of directories
dirnames = []
for filepath in filepaths:
    dirname=filepath.split(os.sep)[-2]
    dirnames.append(dirname)


# x from first file
x,y=np.array(read_xvg(filepaths[0]))

# y from all files
y=[]
for filepath in filepaths:
    xp,yp = read_xvg(filepath)
    y.append(yp)
y=np.array(y)

# save xy
xy=np.append(np.array([x]),y,axis=0)
xy=xy.transpose()


with open(outfile, 'w') as writer:
    #FIXME: The header is hard codec.
    # write header
    writer.write('@    title "Cumulative total radial average charge"' + '\n')
    writer.write('@    xaxis  label "Distance [nm]"' + '\n')
    writer.write('@    yaxis  label "Charge [C]"' + '\n')
    writer.write('@TYPE xy' + '\n')

    # write legends
    for n in range(len(dirnames)):
        writer.write('@ s' + str(n) + ' legend "' + dirnames[n] + '"\n')

    # write xy
    for row in xy:
        line='%12.4f' % row[0]
        for v in row[1:]:
            line = line + ' %12e' % v
        writer.write(line + '\n')


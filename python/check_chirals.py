#!/usr/bin/env python2

import argparse 
import MDAnalysis


### defaults
resname="POPC"


### Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-c", help='coordinate file', required=True)
parser.add_argument("-p", help='topology file')
args = parser.parse_args()

if args.c:
      coordinatefile = args.c
if args.p:
      topologyfile = args.p



### Read atom masses from .itp
topologystream = open(topologyfile, 'r')
mass = dict()

# skip stuff  before [ atoms ] section
while topologystream.readline().find('[ atoms ]') < 0:
    continue

# skip header
topologystream.readline()

# write [ atoms ] section to tempitp
while True:

    line = topologystream.readline()

    # stop after [ atoms ] section
    if line.startswith('['):
        break
    
    # skip comments
    if line.startswith(';'):
        continue

    # parse line
    columns = line.split()
    if len(columns) < 5:
        continue

    atomname = columns[4]
    atommass = columns[7]
    mass[atomname] = round(atommass)
topologystream.close()







#u = MDAnalysis.Universe(coordinatefile) 


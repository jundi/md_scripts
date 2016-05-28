#!/usr/bin/python

import sys
import os
import argparse
from argparse import RawTextHelpFormatter

#--------------------------
# default input/output files
oldmapfile = 'chol.charmm36.map'
newmapfile = 'chol.opls.map'
itpfile = 'chol-opls.itp'
transfile = 'chol-charmm2opls.trans'



#--------------------------------
# commandline parser
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description=
"""This script translates backward mapping files (see citation below) for
different forcefields. Files you need are:

1. Mapping file for some forcefield
2. The .itp file of the forcefield you want to use
3. A translation file which is just a plain text file that has two colums.
   First column should have the atom names of the FF you want to translate and
   second column should have the corresponding atom names of FF you want to
   use.

Going Backward: A Flexible Geometric Approach to Reverse Transformation from
Coarse Grained to Atomistic Models,
Tsjerk A. Wassenaar, Kristyna Pluhackova, Rainer A. Böckmann, Siewert J.
Marrink, and D. Peter Tieleman,
Journal of Chemical Theory and Computation 2014 10 (2), 676-690"""
)
parser.add_argument("-i", help='old map file')
parser.add_argument("-p", help='new topology file')
parser.add_argument("-o", help='output file')
parser.add_argument("-t", help='translation file')
args = parser.parse_args()

if args.i:
  oldmapfile = args.i
if args.p:
  itpfile = args.p
if args.o:
  newmapfile = args.o
if args.t:
  transfile = args.t

if os.path.exists(newmapfile):
  print(newmapfile + ' already exists.')
  exit(0)



#--------------------------
# load translation
transstream = open(transfile, 'r')
newtoold = {}
oldtonew = {}

for line in transstream:

  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 2:
    continue

  newtoold[columns[1]] = columns[0]
  oldtonew[columns[0]] = columns[1]


transstream.close()



#--------------------------
# load atom order from itp
itpstream = open(itpfile, 'r')
order = []

# skip stuf  before [ atoms ] section
while itpstream.readline().find('[ atoms ]') < 0:
  continue

# skip header
itpstream.readline()

# write [ atoms ] section to tempitp
while True:
  line = itpstream.readline()
  if line.startswith('['):
    break
  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 5:
    continue

  atom = columns[4]
  order.append(atom)

itpstream.close()



#--------------------------
# load cg-beads from old map file
oldmapstream = open(oldmapfile, 'r')
beads = {}

# skip stuf  before [ atoms ] section
while oldmapstream.readline().find('[ atoms ]') < 0:
  continue

# save beads
while True:
  line = oldmapstream.readline()
  if line.startswith('['):
    break
  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 3:
    continue

  atom = columns[1]
  beadset = columns[2:]
  beads[atom] = beadset

oldmapstream.close()




#--------------------------
# write new file
oldmapstream = open(oldmapfile, 'r')
newmapstream = open(newmapfile, 'w')

# Until [ atoms ] section, just copy everything
while True:
  line = oldmapstream.readline()
  newmapstream.write(line)
  if line.find('atoms') > 0:
    break
  if len(line) < 1:
    print('[ atoms ] section not found.')
    sys.exit(0)


# write [ atoms ] section
num = 1
for atom in order:
   
  if atom in newtoold:
    oldatom = newtoold[atom]
    beadstr = ' '.join(beads[oldatom])
  else:
    beadstr = ''

  atomstr = "{0:5}{1:6}".format(str(num), atom)
  newline = atomstr+ '   ' + beadstr + '\n'
  newmapstream.write(newline)
  num = num + 1


# find next section after [ atoms ]
while True:
  line = oldmapstream.readline()
  if line.startswith(';') or line.startswith('[') or line.isspace():
    newmapstream.write(line)
    break

# Switch atom names in the last sections
while True:

  line = oldmapstream.readline()

  if line.startswith(';') or line.startswith('[') or line.isspace():
    newmapstream.write(line)
    continue

  if len(line) < 1:
    print('end of file')
    break

  newline=''
  oldatoms = line.split()
  for a in oldatoms:
    newline = newline + oldtonew[a] + ' '

  newmapstream.write(newline + '\n')


oldmapstream.close()
newmapstream.close()

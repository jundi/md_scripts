#!/usr/bin/python
#
# An exercise script written in a python course in summer 2014.
# A smart way to modify .gro file with python.
# 

import sys 		# sys.argv
import MDAnalysis

file_in = sys.argv[1] 
file_out = 'modified.gro'
coord_trans = [1,2,3] 
box_scale = [2,2,2,1,1,1]	# dimensions and angles
selection = 'not resname SOL'

# load .gro
u = MDAnalysis.Universe(file_in) 

# select non waters
s = u.selectAtoms(selection)

# translate coordinates
s.translate(coord_trans)

# double the box size
#u.dimensions = u.dimensions * box_scale

# write new file
s.write(file_out)

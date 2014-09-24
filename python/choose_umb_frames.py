#!/usr/bin/python

import argparse
import os
import numpy
import shutil



# defaults
#--------------------------------
traj = "../pull/run/traj.trr"
index = "index.ndx"
struct = "../pull/run/topol.tpr"
top = "../topology/topol.top"
mdp = "umbrella.mdp"
window = 0.1
do_g_dist = 1



# commandline parser
#--------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-s", help=struct)
parser.add_argument("-f", help=traj)
parser.add_argument("-n", help=index)
parser.add_argument("-p", help=top)
parser.add_argument("--skip_g_dist", help='Dont calculate distances.', action='store_true')
args = parser.parse_args()

if args.s:
  struct = args.s
if args.f:
  traj = args.f
if args.n:
  index = args.n
if args.p:
  top = args.p
if args.skip_g_dist:
  do_g_dist = 0

# absolute paths:
traj=os.path.abspath(traj)
struct=os.path.abspath(struct)
index=os.path.abspath(index)
mdp=os.path.abspath(mdp)
top=os.path.abspath(top)



# g_dist
#--------------------------------
if do_g_dist > 0:
  g_dist_command = ('g_dist' 
      + ' -f ' + traj
      + ' -s ' + struct 
      + ' -n ' + index
      )

  os.system(g_dist_command)



# find frames
#--------------------------------
frame_list = ([])
dists = numpy.loadtxt('dist.xvg',comments='@',skiprows=8,usecols=(0,4))

# Maximum distance
maxdist = dists[:,1].max()

itr = 0
while itr < maxdist:
  diff = abs( dists[:,1] - itr )
  indice = diff.argmin()
  time = dists[indice,0]
  dist = dists[indice,1]
  frame_list.append([time,dist])
  print str(itr) + '   ' + str(dist)
  itr = itr + window


# extract .gro-files and create run input files
#--------------------------------
for itr in frame_list:

  time = int(round(itr[0]))
  dist = round(itr[1],4)
  foldername = str(dist)
  filename = 'confin.gro'
  gromplog = 'grompp.log'

  # backup old directory 
  backup_num = 0
  backup_folder = foldername
  while os.path.exists(backup_folder):
    print(backup_folder + ' exists')
    backup_num = backup_num + 1
    backup_folder = '#' + foldername + '.' + str(backup_num) + '#'
  
  if backup_num > 0:
    shutil.move(foldername, backup_folder)


  # create new directory
  os.mkdir(foldername)
  os.chdir(foldername)


  # dump gro with trjconv
  trjconv_command = ('trjconv' 
      + ' -dump ' + str(time)
      + ' -o ' + filename
      + ' -b ' + str( max( [time - 1, 0] ) )
      + ' -f ' + traj
      + ' -s ' + struct 
      )
  os.system('echo 0 | ' + trjconv_command)


  # grompp
  grompp_command = ('grompp' 
      + ' -c ' + filename
      + ' -p ' + top
      + ' -f ' + mdp
      + ' -n ' + index
      + ' 2> grompp.log'
      )
  os.system(grompp_command)


  # back to main directory
  os.chdir('../')


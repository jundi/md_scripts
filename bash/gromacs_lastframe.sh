#!/bin/bash
#
# Prints last timestamp of xtc/trr/cpt file
#
# Usage:
#
# bash gromacs_lastframe.sh trajectory.xtc
#

cptfile=$1
lastframe=$(gmx check -f $cptfile 2>&1 | grep "Last frame" | awk '{print $NF}')
if [[ -z  $lastframe ]]; then
  lastframe=$(gmx check -f $cptfile 2>&1 | grep "Reading frame" | tail -n1 | awk '{print $NF}')
fi
echo $lastframe

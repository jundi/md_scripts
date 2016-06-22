#!/bin/bash
#
# Looks for GROMACS state.cpt files and prints the time steps.
#

for state in $(find -name "state.cpt"); do
  lframe=$(gmx check -f $state 2>&1 | grep "Last frame")
  tstep=$(echo $lframe | awk '{print $5}')
  folder=$(dirname $state)
  echo "$folder   $tstep"
done | sort

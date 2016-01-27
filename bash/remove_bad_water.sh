#!/bin/bash
# 
# Removes water molecules around selected atom group.
#

### defaults
input="conf.gro"
output="out.gro"
structure=""
index=""
new_index="badwater.ndx"
radius="0.5"
group="Protein"

### arguments
while [[ $# -gt 0 ]]; do    
  case "$1" in
    -f)
      input=$(readlink -f $2)
      shift
      ;;
    -o)
      output=$(readlink -f $2)
      shift
      ;;
    -s)
      structure=$(readlink -f $2)
      shift
      ;;
    -n)
      index="-n $(readlink -f $2)"
      shift
      ;;
    -on)
      new_index=$(readlink -f $2)
      shift
      ;;
    -g)
      group="$2"
      shift
      ;;
  esac
  shift       
done

if [[ -z "$structure" ]]; then
  structure=$input
fi


### create selection
selection="not (group Water and (same resid as (within $radius of group $group)))"
g_select -on $new_index -f $input -s $structure $index -select "$selection"


### write new file
editconf -n $new_index -f $input -o $output



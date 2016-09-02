#!/bin/bash
set -e


help="
Use gmx select (GROMACS 5.1.2) Create index.ndx file for bilayer leaflets. The original groups of index.ndx file are split based on z-coordinate of center of bilayer. Molecules whose COM z-coordinate is smaller than bilayer center z-coordinate belong to leaflet A. Other molecules belong to leaflet B. Groups of two leaflets are saved in file leaflet_A.ndx and leaflet_B.ndx.\n
\n
Example:\n
\t$(basename $0) -n index.ndx -s topol.tpr -b 3.982 \n
Options:\n
\t-n \t Original index file \n
\t-s \t Topology file \n
\t-b \t z-coordinate of membrane center (nm) \n
"

while [[ $# -gt 0 ]]; do    
  case "$1" in
    -s)
      topol=$(readlink -f $2)
      shift
      ;;
    -n)
      index=$(readlink -f $2)
      shift
      ;;
    -b)
      boxcenter=$2
      shift
      ;;
    *)
      echo -e $help
      exit 2
      ;;
  esac
  shift       
done


# check args
if [[ ! -e ${topol} || ! -e ${index} || -z ${boxcenter} ]]; then
  echo -e ${help}
  exit 3
fi


# set operators for leaflet A and leaflet B
declare -A operator
operator[A]="<"
operator[B]=">"


# get all groups from old index file
groups=$(grep "\[" ${index} | awk '{print $2}')


# create .ndx for leaflets A and B
for leaflet in A B; do

  # create selection string
  selections=""
  for group in ${groups}; do
    selections="${selections} group \"${group}\" and z ${operator[${leaflet}]} ${boxcenter};"
  done

  # create new index file
  gmx select -n ${index} -s ${topol} -selrpos mol_com -select "${selections}" -on leaflet_${leaflet}.ndx

  # rename groups in new index file
  sed -i "s/\"_and_z_${operator[${leaflet}]}_${boxcenter}_f0_t0.000 \]/ \]/" leaflet_${leaflet}.ndx
  sed -i "s/\[ group_\"/\[ /" leaflet_${leaflet}.ndx

done

#!/bin/bash

####################A######################
#
# Input:
# -n index.ndx
# -s topol.tpr
# -f traj.xtc
# -b first frame to use
# 
# analysis to be done
# * rdf
# * sorient
# * order
# * rms
# * potential
# * mindist (very slow)
# * sas (very slow)
# 
##########################################

####################A######################
#
# Index file should contain these groups: 
# * System
# * CO (Cholesteryl oleate)
# * POPC
# * Protein
# * NA
# * CL
# * Water
# * Lipids (All lipids)
# * HDL (Lipids and Proteins)
# * POPC_P (POPC Phosphate P-)
# * POPC_N (POPC Choline N+)
# * POPC_Protein (POPC and Protein)
#
##########################################


############
# defaults #
############
traj="traj.xtc"
structure="topol.tpr"
index="index.ndx"
begin=0
tasks=() # empty array


####################
# input parameters #
####################
if [[ $# -lt 1 ]]; then
  echo "Not enough input parameteres."
  exit 1
fi
while [[ $# -gt 0 ]]; do    
  case "$1" in
    -s)
      structure="$2"
      shift
      ;;
    -f)
      traj="$2"
      shift
      ;;
    -n)
      index="$2"
      shift
      ;;
    -b)
      begin="$2"
      shift
      ;;
    rdf|sorient|order|rms|potential|mindist|sas)
      tasks+=("$1")
      ;;
    *)
      echo "RTFM"
      exit 2
      ;;
  esac
  shift       
done



########
# main #
########
main() {

  for task in ${tasks[@]}
  do
    echo -e "Calculating $task..."
    $task  >"$task"".log" 2> "$task""2.log" && echo -e "$task done!"
  done

  wait
  echo -e "All tasks completed."

}



#######
# RDF #
#######
rdf() {

  workdir=g_rdf
  mkdir -p $workdir
  cd $workdir

  ref_group="Lipids"
  groups="CO POPC Protein NA CL Water POPC_P POPC_N"
  ng=8

  echo "$ref_group $groups" | g_rdf -f ../$traj -n ../$index -s ../$structure  -b $begin -rdf atom -com -ng 8

  cd ..
}



########
# RMSD #
########
rms() {

  workdir=g_rms
  mkdir -p $workdir
  cd $workdir

  ref_group="Lipids"
  groups="CO POPC Protein Lipids HDL"

  echo "$ref_group $groups" | g_rms -f ../$traj -n ../$index -s ../$structure -ng 5 -what rmsd 

  cd ..
}



###################
# ORDER PARAMETER #
###################


order() {

  workdir=g_order
  mkdir -p $workdir
  cd $workdir

  # create selection string for g_select
  resname="POPC"
  # For some reason one extra group has to be added to the index file (when
  # using radial-option).  Otherwise order parameter of last atom will be
  # missing
  palmitoyl=(C36 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C52) 
  select=''
  for atom in ${palmitoyl[@]}; do
    select="$select name $atom and resname $resname;"
  done

  # create index file for palmitoyl
  ndx_tail="$resname-palmitoyl.ndx"
  g_select -s ../$structure -select "$select" -on $ndx_tail

  # Reference group
  ref_group="Lipids"
  # Output file
  order_xvg="order-$resname""_palmitoyl.xvg"

  echo "$ref_group" | g_order -f ../$traj -nr ../$index -s ../$structure  -b $begin -n $ndx_tail -radial -permolecule -o $order_xvg

  cd ..
}



###########################
# ELECTROSTATIC POTENTIAL #
###########################
potential() {

  workdir=g_H_potential
  mkdir -p $workdir
  cd $workdir

  # I have a different version of GROMACS for using g_H_potential
  binsize=0.05
  # Reference group
  ref_group="Lipids"
  groups=(CO POPC Protein NA CL Water System POPC_N POPC_P)

  for group in ${groups[@]}
  do
    mkdir $group
    cd $group
    select="com of group $ref_group pbc; group $group"

    g_H_potential -f ../../$traj -n ../../$index -s ../../$structure  -b $begin -geo Radial -bin_size $binsize -select "$select" &
    cd ..

    while [[ $(jobs | wc -l) -gt 6 ]]; do
      sleep 5
    done

  done

  cd ..
}


#####################
# WATER ORIENTATION #
#####################
sorient() {

  workdir=g_sorient
  mkdir -p $workdir
  cd $workdir

  # Reference group
  ref_group="Lipids"
  # water group
  group="Water"

  # Calculate sorient for 0.1nm slices
  rstep=0.1
  rmin=3
  rmax=$(echo "$rmin+$rstep" | bc -l)
  rmaxmax=7
  dt=1000

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    echo "$ref_group $group" | g_sorient -f ../../$traj -n ../../$index -s ../../$structure -b $begin -com -rmin $rmin -rmax $rmax -dt $dt &
    cd ..

    while [[ $(jobs | wc -l) -gt 6 ]]; do
      sleep 5
    done

    #next slice:
    rmin=$rmax
    rmax=$(echo "$rmin+$rstep" | bc -l)

  done

  cd ..

}




############
# CONTACTS #
############
mindist() {

  workdir=g_mindist
  mkdir -p $workdir
  cd $workdir

  dist=0.25
  ref_groups=(CO POPC Protein Lipids HDL POPC_Protein)
  groups="NA CL Water"
  dt=1000 # 1ns

  for ref_group in ${ref_groups[@]}; do
    echo "$ref_group $groups" | g_mindist -f ../$traj -n ../$index -s ../$structure -group -ng 3 -dt $dt -od "$ref_group-mindist.xvg" -on "$ref_group-numcount.xvg" -d $dist &

    while [[ $(jobs | wc -l) -gt 6 ]]; do
      sleep 5
    done

  done

  cd ..
}




########
# SASA #
########
sas() {

  workdir=g_sas
  mkdir -p $workdir
  cd $workdir

  ref_group="HDL"
  groups=(CO POPC Protein Lipids HDL)
  dt=1000 # 1ns

  for group in  ${groups[@]}; do
    echo "$ref_group $group" | g_sas -f ../$traj -n ../$index -s ../$structure -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb -dt $dt &
  done

  while [[ $(jobs | wc -l) -gt 6 ]]; do
    sleep 5
  done

  cd ..
}



#####################
# run main function #
#####################
main

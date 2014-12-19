#!/bin/bash

####################A######################
# Input:
# -n index.ndx
# -s topol.tpr
# -f traj.xtc
# -b first frame to use (ps)
# -j max jobs
# 
# tasks:
# * rdf
# * sorient
# * order
# * rms
# * potential
# * mindist (very slow)
# * sas (very slow)
#
# example:
# $ sh batch-analyse.sh -n index.ndx -s topol.tpr -f traj.xtc -b 100000 -j 4 rdf sorient sas
#
##########################################

####################A######################
# Index file should contain these groups: 
#
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
maxjobs=2 # max parallel jobs

#############
# variables #
#############
tasks=() # array to store tasks
jobid=$BASHPID # jobid for jobcontrol


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
    -j)
      maxjobs="$2"
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
    $task  >"$task"".log" 2> "$task""2.log"
  done

  wait
  echo -e "All tasks completed."

}



###############
# job limiter #
###############
waitjobs() {

  jobnum=$(ps -o pgid | grep $jobid | wc -l)
  let jobnum=$jobnum-5

  while [[ $jobnum -ge $maxjobs ]]
  do
    sleep 10

    jobnum=$(ps -o pgid | grep $jobid | wc -l)
    let jobnum=$jobnum-5

  done
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

  # g_rdf
  echo "$ref_group $groups" | g_rdf -f ../$traj -n ../$index -s ../$structure  -b $begin -rdf atom -com -ng 8 &

  # wait until other jobs finish
  waitjobs

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

  # g_rms
  echo "$ref_group $groups" | g_rms -f ../$traj -n ../$index -s ../$structure -ng 5 -what rmsd &

  # wait until other jobs finish
  waitjobs

  cd ..
}



###################
# ORDER PARAMETER #
###################


order() {

  workdir=g_order
  mkdir -p $workdir
  cd $workdir



  ### create index file for palmitoyl
  resname="POPC"
  # For some reason one dummy group has to be added to the index file (when
  # using radial-option).  Otherwise order parameter of last atom will be
  # missing
  palmitoyl=(C36 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C52) 
  select=''

  for atom in ${palmitoyl[@]}; do
    select="$select name $atom and resname $resname;"
  done

  ndx_tail="$resname-palmitoyl.ndx"
  g_select -s ../$structure -select "$select" -on $ndx_tail



  # Reference group
  ref_group="Lipids"
  # Output file
  order_xvg="order-$resname""_palmitoyl.xvg"

  # g_order
  echo "$ref_group" | g_order -f ../$traj -nr ../$index -s ../$structure  -b $begin -n $ndx_tail -radial -permolecule -o $order_xvg &

  # wait until other jobs finish
  waitjobs

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

    # g_H_potential
    g_H_potential -f ../../$traj -n ../../$index -s ../../$structure  -b $begin -geo Radial -bin_size $binsize -select "$select" &

    # wait until other jobs finish
    waitjobs
    cd ..

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
  dt=1000 # 1ns

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    # g_sorient
    echo "$ref_group $group" | g_sorient -f ../../$traj -n ../../$index -s ../../$structure -b $begin -com -rmin $rmin -rmax $rmax -dt $dt &

    # wait until other jobs finish
    waitjobs


    #next slice:
    rmin=$rmax
    rmax=$(echo "$rmin+$rstep" | bc -l)
    cd ..

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
  ref_groups=(HDL Protein Lipids CO POPC POPC_Protein)
  groups="NA CL Water"
  dt=1000 # 1ns

  for ref_group in ${ref_groups[@]}; do

    # g_mindist
    echo "$ref_group $groups" | g_mindist -f ../$traj -n ../$index -s ../$structure -group -ng 3 -dt $dt -od "$ref_group-mindist.xvg" -on "$ref_group-numcount.xvg" -d $dist &

    # wait until other jobs finish
    waitjobs

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

    # g_sas
    echo "$ref_group $group" | g_sas -f ../$traj -n ../$index -s ../$structure -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb -dt $dt & 

    # wait until other jobs finish
    waitjobs

  done

  cd ..
}



#####################
# run main function #
#####################
main

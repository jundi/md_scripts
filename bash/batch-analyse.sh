#!/bin/bash

####################A######################
# Input:
# -n index.ndx
# -nn index_nowater.ndx
# -s topol.tpr
# -sn topol_nowater.tpr
# -f traj.xtc
# -fn traj_nowater.xtc
# -b first frame to use (ps)
# -dt skip frames
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
# * dssp 
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
# * CHO (Cholesteryl oleate)
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
# input files
traj=$(readlink -f traj.xtc)
structure=$(readlink -f topol.tpr)
index=$(readlink -f index.ndx)
# input files without water (used when water molecules are not needed for
# analysis)
traj_nw=""
structure_nw=""
index_nw=""
# other parameters
begin=0   # first timestep to be used
dt=-1     # skip frames
maxjobs=2 # max parallel jobs


####################
# global variables #
####################
tasks=()       # array to store tasks
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
      structure=$(readlink -f $2)
      shift
      ;;
    -sn)
      structure_nw=$(readlink -f $2)
      shift
      ;;
    -f)
      traj=$(readlink -f $2)
      shift
      ;;
    -fn)
      traj_nw=$(readlink -f $2)
      shift
      ;;
    -n)
      index=$(readlink -f $2)
      shift
      ;;
    -nn)
      index_nw=$(readlink -f $2)
      shift
      ;;
    -b)
      begin="$2"
      shift
      ;;
    -dt)
      dt="$2"
      shift
      ;;
    -j)
      maxjobs="$2"
      shift
      ;;
    rdf|sorient|order|rms|potential|mindist|sas|dssp)
      tasks+=("$1")
      ;;
    *)
      echo "RTFM"
      exit 2
      ;;
  esac
  shift       
done

# if non-water input files are not given, use input files with water
if [[ -z $traj_nw ]]; then
  traj_nw=$traj
fi
if [[ -z $structure_nw ]]; then

  structure_nw=$structure
fi
if [[ -z $index_nw ]]; then
  index_nw=$index
fi



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
  groups="CO CHO POPC Protein NA CL Water POPC_P POPC_N"
  ng=9

  # g_rdf
  echo "$ref_group $groups" | g_rdf -f $traj -n $index -s $structure  -b $begin -rdf atom -com -ng $ng  -dt $dt -cn & 

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
  groups="CO CHO POPC Protein Lipids HDL"

  # g_rms
  echo "$ref_group $groups" | g_rms -f $traj_nw -n $index_nw -s $structure_nw -ng 6 -what rmsd  -dt $dt & 

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
  palmitoyl_carbons=(C36 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C52) 
  palmitoyl_select=''
  for atom in ${palmitoyl_carbons[@]}; do
    palmitoyl_select="$palmitoyl_select name $atom and resname $resname;"
  done
  palmitoyl_ndx="${resname}_palmitoyl.ndx"
  g_select -s $structure_nw -select "$palmitoyl_select" -on $palmitoyl_ndx

  # oleyl chain
  oleyl_carbons=(C15 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C32 C33 C33) 
  oleyl_select=''
  for atom in ${oleyl_carbons[@]}; do
    oleyl_select="$oleyl_select name $atom and resname $resname;"
  done
  oleyl_ndx="${resname}_oleyl.ndx"
  g_select -s $structure_nw -select "$oleyl_select" -on $oleyl_ndx

  # unsaturated part of oleyl chain
  oleyl_unsat_carbons=(C23 C24 C25 C26 C26) 
  oleyl_unsat_select=''
  for atom in ${oleyl_unsat_carbons[@]}; do
    oleyl_unsat_select="$oleyl_unsat_select name $atom and resname $resname;"
  done
  oleyl_unsat_ndx="${resname}_oleyl_unsat.ndx"
  g_select -s $structure_nw -select "$oleyl_unsat_select" -on $oleyl_unsat_ndx



  # Reference group
  ref_group="Lipids"

  # Output files
  order_palmitoyl="order_${resname}_palmitoyl.xvg"
  sliced_palmitoyl="sliced_${resname}_palmitoyl.xvg"
  order_oleyl="order_${resname}_oleyl.xvg"
  sliced_oleyl="sliced_${resname}_oleyl.xvg"
  order_oleyl_unsat="order_${resname}_oleyl_unsat.xvg"
  sliced_oleyl_unsat="sliced_${resname}_oleyl_unsat.xvg"

  # g_order
  echo "$ref_group" | g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $begin -n $palmitoyl_ndx -radial -permolecule -o $order_palmitoyl -os $sliced_palmitoyl -dt $dt &
  # wait until other jobs finish
  waitjobs
  echo "$ref_group" | g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $begin -n $oleyl_ndx -radial -permolecule -o $order_oleyl -os $sliced_oleyl -dt $dt &
  # wait until other jobs finish
  waitjobs
  echo "$ref_group" | g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $begin -n $oleyl_unsat_ndx -radial -permolecule -o $order_oleyl_unsat -os $sliced_oleyl_unsat -dt $dt &
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
  groups_nw=(CO CHO POPC Protein NA CL POPC_N POPC_P)
  groups_water=(Water System)

  for group in ${groups_nw[@]}
  do
    mkdir $group
    cd $group
    select="com of group $ref_group pbc; group $group"

    # g_H_potential
    g_H_potential -f $traj_nw -n $index_nw -s $structure_nw  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt &

    # wait until other jobs finish
    waitjobs
    cd ..

  done

  for group in ${groups_water[@]}
  do
    mkdir $group
    cd $group
    select="com of group $ref_group pbc; group $group"

    # g_H_potential
    g_H_potential -f $traj -n $index -s $structure  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt &

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

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    # g_sorient
    echo "$ref_group $group" | g_sorient -f $traj -n $index -s $structure -b $begin -com -rmin $rmin -rmax $rmax -dt $dt &

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

  dist=0.35
  ref_groups=(HDL Protein Lipids CO CHO POPC POPC_Protein)


  groups="NA CL"
  for ref_group in ${ref_groups[@]}; do

    # g_mindist
    echo "$ref_group $groups" | g_mindist -f $traj_nw -n $index_nw -s $structure_nw -group -ng 3 -dt $dt -od "$ref_group-mindist.xvg" -on "$ref_group-numcount.xvg" -d $dist &

    # wait until other jobs finish
    waitjobs

  done

  groups="Water"
  for ref_group in ${ref_groups[@]}; do

    # g_mindist
    echo "$ref_group $groups" | g_mindist -f $traj -n $index -s $structure -group -ng 3 -dt $dt -od "$ref_group-mindist.xvg" -on "$ref_group-numcount.xvg" -d $dist &

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
  groups=(CO CHO POPC Protein Lipids HDL)

  for group in  ${groups[@]}; do

    # g_sas
    echo "$ref_group $group" | g_sas -f $traj_nw -n $index_nw -s $structure_nw -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb -dt $dt & 

    # wait until other jobs finish
    waitjobs

  done

  cd ..
}



########
# DSSP #
########
dssp() {

  workdir=do_dssp
  mkdir -p $workdir
  cd $workdir

  chains=(Chain_1 Chain_2 Chain_3 Chain_4)
  for chain in ${chains[@]}; do

    mkdir -p $chain
    cd $chain

    # do_dssp
    echo "$chain" | do_dssp -f $traj_nw -n $index_nw -s $structure_nw -dt $dt & 

    # wait until other jobs finish
    waitjobs

    cd ..

  done

  cd ..
}



#####################
# run main function #
#####################
main

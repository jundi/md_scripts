#!/bin/bash

####################A######################
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
##########################################

############
# defaults #
############
#traj="traj.xtc"
structure="topol.tpr"
index="index.ndx"
begin=0

####################
# input parameters #
####################
if [[ $# -lt 1 ]]; then
  echo "not enough input parameteres"
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
    *)
      analys="$1"
  esac
  shift       
done



########
# main #
########
main() {
  # run one of the functions:
  #rdf
  #rms
  #potential
  #sorient
  #sas
  #mindist
  $analys
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

  echo "$ref_group $groups" | g_rdf -f ../$traj -n ../$index -s ../$structure  -b $begin -rdf atom -com -ng 7

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

  # which are the C-atoms of the tail?
  first_atom=38
  last_atom=52
  resname="POPC"

  # create index file containing the tails using a script
  /home/mikkolah/koodi/bash/g_select_tails.sh "$first_atom" "$last_atom" "$resname" "../$structure"

  # Reference group
  ref_group="Lipids"

  echo "$ref_group" | g_order -f ../$traj -nr ../$index -s ../$structure  -b $begin -n $resname-$first_atom-$last_atom.ndx -radial -permolecule -o $resname-$first_atom-$last_atom.xvg

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
    g_H_potential -f ../../$traj -n ../../$index -s ../../$structure  -b $begin -geo Radial -bin_size $binsize -select "$select"
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

  # Calculate sorient for 0.5nm slices
  rstep=0.5
  rmin=0
  rmax=$rstep
  rmaxmax=10

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    echo "$ref_group $group" | g_sorient -f ../../$traj -n ../../$index -s ../../$structure -b $begin -com -rmin $rmin -rmax $rmax
    cd ..

    #next slice:
    rmin=$rmax
    rmax=$(echo "$rmax+$rstep" | bc -l)

  done

  cd ..

}




###########
# CONTATS #
###########
mindist() {

  workdir=g_mindist
  mkdir -p $workdir
  cd $workdir

  dist=0.25
  ref_groups=(CO POPC Protein Lipids HDL POPC_Protein)
  groups="NA CL Water"
  dt=1000 # 1ns

  for ref_group in ${ref_groups[@]}; do
    echo "$ref_group $groups" | g_mindist -f ../$traj -n ../$index -s ../$structure -group -ng 3 -dt $dt -od "$ref_group-mindist.xvg" -on "$ref_group-numcount.xvg" -d $dist
  done

  cd ..
}




########
# SASA #
########
sas() {

  echo "jee"
  workdir=g_sas
  mkdir -p $workdir
  cd $workdir

  ref_group="HDL"
  groups=(CO POPC Protein Lipids HDL)

  for group in  ${groups[@]}; do
    echo "$ref_group $group" | g_sas -f ../$traj -n ../$index -s ../$structure -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb
  done

  cd ..
}



#####################
# run main function #
#####################
main

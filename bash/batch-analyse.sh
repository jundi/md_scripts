#!/bin/bash

usage="\n
Usage: \n
\t$(basename $0) [OPTION...] [TASK1,TASK2,...] \n
\n
Example: \n
\t$ sh batch-analyse.sh -n index.ndx -s topol.tpr -f traj.xtc -b 100000 -j 4 rdf sorient sas \n
\n
Options: \n
\t-n \t index.ndx \n
\t-nn \t index_nowater.ndx \n
\t-s \t topol.tpr \n
\t-sn \t topol_nowater.tpr \n
\t-f \t traj.xtc \n
\t-fn \t traj_nowater.xtc \n
\t-b \t first frame to use (ps) \n
\t-dt \t skip frames \n
\t-j \t max parallel jobs \n
\n
Tasks: \n
\tdssp \n
\tgyrate \n
\torder \n
\tmindist (water needed)  \n
\tpotential (water needed) \n
\trdf \n
\trms \n
\tsas \n
\tsorient (water needed) \n
\n
Index file should contain these groups: \n
\tSystem \n
\tCO (Cholesteryl oleate) \n
\tCHO (Cholesterol) \n
\tPOPC and/or DPPC\n
\tProtein \n
\tNA \n
\tCL \n
\tWater \n
\tMonolayer (CHO and POPCs) \n
\tLipids (All lipids) \n
\tHDL (Lipids and Proteins) \n
\tPOPC_P and/or DPPC_P (PL Phosphate P-) \n
\tPOPC_N and/or DPPC_N (PL Choline N+) \n
\tPOPC_Protein and/or DPPC_Protein (PL and Protein) \n
\tChain_1 \n
\tChain_2 \n
\tChain_3 \n
\tChain_4 \n
"



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
    -h)
      echo -e $usage
      exit 0
      ;;
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
    rdf|sorient|order|rms|potential|mindist|sas|dssp|gyrate)
      tasks+=("$1")
      ;;
    *)
      echo -e $usage
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





#######
# RDF #
#######
rdf() {

  workdir=g_rdf
  bin=0.02
  mkdir -p $workdir
  cd $workdir

  ref_group="Lipids"
  groups=""
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "Water" "POPC_P" "POPC_N" "DPPC_P" "DPPC_N")
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups="${groups} ${g}"
    fi
  done
  ng=$(echo $groups | wc -w)

  # rdf
  echo "$ref_group $groups" | sem -j $maxjobs g_rdf -f $traj -n $index -s $structure  -b $begin -bin $bin -rdf atom -com -ng $ng  -dt $dt -cn rdf_cn.xvg -o rdf.xvg 

  # nonorm
  echo "$ref_group $groups" | g_rdf -f $traj -n $index -s $structure  -b $begin -bin $bin -rdf atom -com -ng $ng  -dt $dt -cn nonorm_cn.xvg -o nonorm.xvg 

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
  groups=""
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "Monolayer" "Lipids" "HDL")
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups="${groups} ${g}"
    fi
  done
  ng=$(echo $groups | wc -w)

  # g_rms
  echo "$ref_group $groups" | g_rms -f $traj_nw -n $index_nw -s $structure_nw -ng $ng -what rmsd  -dt $dt 

  cd ..
}



###################
# ORDER PARAMETER #
###################


order() {

  workdir=g_order
  mkdir -p $workdir
  cd $workdir


  tailnames=("POPC_SN1" "POPC_SN2" "POPC_SN2_unsat" "DPPC_SN1" "DPPC_SN2")
  for tn in ${tailnames[@]}; do

    # TAIL ATOMS
    # For some reason (bug?) one dummy group has to be added to the index file
    # (when using radial-option).  Otherwise order parameter of last atom will
    # be missing.
    case $tn in
      "POPC_SN1")
	tail=(C36 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C52)
	resname=POPC
	unsat="-nounsat"
	;;
      "POPC_SN2")
	tail=(C15 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C32 C33 C33)
	resname=POPC
	unsat="-nounsat"
	;;
      "POPC_SN2_unsat")
	tail=(C23 C24 C25 C26 C26)
	resname=POPC
	unsat="-unsat"
	;;
      "DPPC_SN1")
	tail=(C34 C36 C37 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C50) 
	resname=DPPC
	unsat="-nounsat"
	;;
      "DPPC_SN2")
	tail=(C15 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C31) 
	resname=DPPC
	unsat="-nounsat"
	;;
      *)
	echo "ERROR: Unknown tail"
	continue
	;;
    esac

    if [[ ! $(grep " $resname " $index) ]]; then
      continue
    fi


    # BUILD INDEX FILE
    select=""
    for atom in ${tail[@]}; do
      select="$select name $atom and resname $resname;"
    done

    tail_ndx="${tn}.ndx"
    g_select -s $structure_nw -select "$select" -on $tail_ndx


    # G_ORDER
    # Reference group
    ref_group="Lipids"
    echo "$ref_group" | g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $begin -n $tail_ndx -radial -permolecule -o "order_$tn" -os "sliced_$tn" -dt $dt $unsat

  done

  cd ..
}



###########################
# ELECTROSTATIC POTENTIAL #
###########################
potential() {

  workdir=g_H_potential
  mkdir -p $workdir
  cd $workdir

  binsize=0.001

  # Reference group
  ref_group="Lipids" 

  # groups which include water
  groups_water=("Water" "System")

  # buil list of groups which do not include water
  groups_nw=()
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "POPC_N" "POPC_P" "DPPC_N" "DPPC_P")
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups_nw+=("$g")
    fi
  done

  # calculations using trajectories WITHOUT water
  for group in ${groups_nw[@]}
  do
    mkdir -p $group
    cd $group
    select="\"com of group $ref_group pbc; group $group\""

    # g_H_potential
    sem -j $maxjobs g_H_potential -f $traj_nw -n $index_nw -s $structure_nw  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt 

    cd ..

  done

  # calculations using trajectories WITH water
  for group in ${groups_water[@]}
  do
    mkdir $group
    cd $group
    select="\"com of group $ref_group pbc; group $group\""

    # g_H_potential
    g_H_potential -f $traj -n $index -s $structure  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt 

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

  # Calculate sorient for 0.3nm slices
  rstep=0.3
  rmin=3
  rmax=$(echo "$rmin+$rstep" | bc -l)
  rmaxmax=7
  cbin=0.05
  rbin=0.05

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    # g_sorient
    echo "$ref_group $group" | g_sorient -f $traj -n $index -s $structure -b $begin -cbin $cbin -rbin $rbin -com -rmin $rmin -rmax $rmax -dt $dt 

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
  ref_groups=()
  group_list=(HDL Protein Lipids Monolayer CO CHO POPC DPPC POPC_Protein DPPC_Protein)
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      ref_groups+=("$g")
    fi
  done


  groups="NA CL"
  for ref_group in ${ref_groups[@]}; do

    # g_mindist
    echo "$ref_group $groups" | sem -j $maxjobs g_mindist -f $traj_nw -n $index_nw -s $structure_nw -group -ng 2 -dt $dt -od "${ref_group}-ions_mindist.xvg" -on "${ref_group}-ions_numcount.xvg" -d $dist

  done

  groups="Water"
  for ref_group in ${ref_groups[@]}; do

    # g_mindist
    echo "$ref_group $groups" | g_mindist -f $traj -n $index -s $structure -group -ng 1 -dt $dt -od "${ref_group}-water_mindist.xvg" -on "${ref_group}-water_numcount.xvg" -d $dist 

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

  # reference group is whole HDL particle or lipid droplet (if there is no
  # protein at all)
  if [[ $(grep " HDL " $index) ]]; then
    ref_group="HDL"
  else
    ref_group="Lipids"
  fi

  # target groups
  groups=()
  group_list=(CO CHO POPC DPPC Protein Monolayer Lipids HDL)
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups+=("$g")
    fi
  done

  for group in  ${groups[@]}; do


    # g_sas
    echo "$ref_group $group" | g_sas -f $traj_nw -n $index_nw -s $structure_nw -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb -dt $dt 

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
    echo "$chain" | do_dssp -f $traj_nw -n $index_nw -s $structure_nw -dt $dt 


    cd ..

  done

  cd ..
}



###################
# GYRATION RADIUS #
###################
gyrate() {

  workdir=g_gyrate
  mkdir -p $workdir
  cd $workdir

  groups=(CO Lipids HDL Protein)

  for group in ${groups[@]}; do

    echo $group | sem -j $maxjobs g_gyrate -f $traj_nw -n $index_nw -s $structure_nw -dt $dt -o ${group}.xvg

  done

  cd ..
}



#####################
# run main function #
#####################
main

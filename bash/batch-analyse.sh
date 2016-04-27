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
\t-e \t last frame to use (ps) \n
\t-r \t HDL radius \n
\t-dt \t skip frames \n
\t-j \t max parallel jobs \n
\t-jn \t max parallel jobs for analysis without water\n
\n
Tasks: \n
\tdssp \n
\tgyrate \n
\torder \n
\torder_blocks \n
\tmindist (water needed)  \n
\tpotential (water needed) \n
\tpotential_blocks (water needed) \n
\trdf (water needed) \n
\trdf blocks (water needed) \n
\trms \n
\tsas \n
\tsorient (water needed) \n
\tsorient_blocks (water needed) \n
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
traj=""
structure=""
index=""
# input files without water (used when water molecules are not needed for
# analysis)
traj_nw=""
structure_nw=""
index_nw=""
# other parameters
begin=0   # first timestep to be used
dt=-1     # skip frames
maxjobs_nw=4 # max parallel jobs
maxjobs=1 # max parallel jobs


####################
# global variables #
####################
tasks=()       # array to store tasks


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
    -e)
      lastframe="$2"
      shift
      ;;
    -r)
      radius="$2"
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
    -jn)
      maxjobs_nw="$2"
      shift
      ;;
    rdf|sorient|order|rms|potential|mindist|sas|dssp|gyrate|potential_blocks|sorient_blocks|order_blocks|rdf_blocks)
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

  sem --wait
  echo -e "All tasks completed."

}





#######
# RDF #
#######
rdf() {

  workdir=g_rdf
  mkdir -p $workdir
  cd $workdir

  bin=0.02
  ref_group="Lipids"
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "Water" "POPC_P" "POPC_N" "DPPC_P" "DPPC_N")

  # Build group list
  groups=""
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups="${groups} ${g}"
    fi
  done
  ng=$(echo $groups | wc -w) # number of groups

  # rdf
  echo "$ref_group $groups" | sem -j $maxjobs g_rdf -f $traj -n $index -s $structure  -b $begin -bin $bin -rdf atom -com -ng $ng  -dt $dt -cn rdf_cn.xvg -o rdf.xvg 

  # nonorm
  echo "$ref_group $groups" | sem -j $maxjobs g_rdf -f $traj -n $index -s $structure  -b $begin -bin $bin -rdf atom -com -ng $ng  -dt $dt -cn nonorm_cn.xvg -o nonorm.xvg 

  cd ..
}

##############
# RDF BLOCKS #
##############
rdf_blocks() {

  workdir=g_rdf_blocks
  mkdir -p $workdir
  cd $workdir

  blocksize=100000
  ref_group="Lipids"
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "POPC_P" "POPC_N" "DPPC_P" "DPPC_N" "Monolayer" "Lipids" "HDL" "Hydrophobic_AA" "Hydrophile_AA" "POPC_C33" "POPC_C52" "DPPC_C31" "DPPC_C50")
  group_list_water=("Water" "not_CO" "not_Lipids" "Water_and_ions")


  groups=""
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      groups="${groups} ${g}"
    fi
  done

  groups_water=""
  for g in ${group_list_water[@]}; do
    if [[ $(grep " $g " $index) ]]; then
      groups_water="${groups_water} ${g}"
    fi
  done


  for group in $groups;do
    mkdir -p $group
    cd $group

    b=1
    while [[ $b -lt $lastframe ]]; do

      let e=${b}+${blocksize}-1
      mkdir -p ${b}-${e}
      cd ${b}-${e}

      #g_rdf
      echo "$ref_group $group" | sem -j $maxjobs_nw g_rdf -f $traj_nw -n $index_nw -s $structure_nw  -b $b -e $e -rdf atom -com -dt $dt -bin 0.02 -o nonorm.xvg -nonorm 
      echo "$ref_group $group" | sem -j $maxjobs_nw g_rdf -f $traj_nw -n $index_nw -s $structure_nw  -b $b -e $e -rdf atom -com -dt $dt -bin 0.02 -o rdf.xvg 

      let b=${b}+${blocksize}
      cd ..
    done
    cd ..
  done

  for group in $groups_water ;do
    mkdir -p $group
    cd $group

    b=1
    while [[ $b -lt $lastframe ]]; do

      let e=${b}+${blocksize}-1
      mkdir -p ${b}-${e}
      cd ${b}-${e}

      #g_rdf
      echo "$ref_group $group" | sem -j $maxjobs g_rdf -f $traj -n $index -s $structure  -b $b -e $e -rdf atom -com -dt $dt -bin 0.02 -o nonorm.xvg -nonorm 
      echo "$ref_group $group" | sem -j $maxjobs g_rdf -f $traj -n $index -s $structure  -b $b -e $e -rdf atom -com -dt $dt -bin 0.02 -o rdf.xvg 

      let b=${b}+${blocksize}
      cd ..
    done

    cd ..
  done

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
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "Monolayer" "Lipids" "HDL")

  # Build group list
  groups=""
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      groups="${groups} ${g}"
    fi
  done
  ng=$(echo $groups | wc -w) # number of groups

  # g_rms
  echo "$ref_group $groups" | sem -j $maxjobs_nw g_rms -f $traj_nw -n $index_nw -s $structure_nw -ng $ng -what rmsd  -dt $dt 

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
  ref_group="Lipids"

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

    if [[ ! $(grep " $resname " $index_nw) ]]; then
      continue
    fi


    # BUILD INDEX FILE
    select=""
    for atom in ${tail[@]}; do
      select="$select name $atom and resname $resname;"
    done

    tail_ndx="${tn}.ndx"
    g_select -s $structure_nw -select "$select" -on $tail_ndx


    # g_order
    echo "$ref_group" | sem -j $maxjobs_nw g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $begin -n $tail_ndx -radial -permolecule -o "order_$tn" -os "sliced_$tn" -dt $dt $unsat

  done

  cd ..
}


##########################
# ORDER PARAMETER BLOCKS #
##########################
order_blocks() {

  workdir=g_order_blocks
  mkdir -p $workdir
  cd $workdir

  blocksize=100000 # 100 ns
  ref_group="Lipids" # Reference group

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


    if [[ ! $(grep " $resname " $index_nw) ]]; then
      continue
    fi


    # BUILD INDEX FILE
    select=""
    for atom in ${tail[@]}; do
      select="$select name $atom and resname $resname;"
    done
    tail_ndx="${tn}.ndx"
    g_select -s $structure_nw -select "$select" -on $tail_ndx


    b=1
    while [[ $b -lt $lastframe ]]; do
      let e=${b}+${blocksize}-1
      mkdir -p ${b}-${e}
      cd ${b}-${e}

      # g_order
      echo "$ref_group" | sem -j $maxjobs_nw g_order -f $traj_nw -nr $index_nw -s $structure_nw  -b $b -e $e -n ../$tail_ndx -radial -permolecule -o "order_$tn" -os "sliced_$tn" -dt $dt $unsat

      cd ..
      let b=${b}+${blocksize}
    done

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
  ref_group="Lipids" # Reference group
  groups_water=("Water" "System") # groups which include water
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "POPC_N" "POPC_P" "DPPC_N" "DPPC_P")

  # build list of groups which do not include water
  groups_nw=()
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      groups_nw+=("$g")
    fi
  done

  # calculations using trajectories WITHOUT water
  for group in ${groups_nw[@]}
  do
    mkdir -p $group
    cd $group
    select="\"com of group $ref_group pbc; group $group\""
      echo $group

    # g_H_potential
    sem -j $maxjobs_nw g_H_potential -f $traj_nw -n $index_nw -s $structure_nw  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt 

    cd ..

  done

  # calculations using trajectories WITH water
  for group in ${groups_water[@]}
  do
    mkdir -p $group
    cd $group
    select="\"com of group $ref_group pbc; group $group\""

    # g_H_potential
    sem -j $maxjobs g_H_potential -f $traj -n $index -s $structure  -b $begin -geo Radial -bin_size $binsize -select "$select" -dt $dt 

    cd ..

  done

  cd ..
}



##################################
# ELECTROSTATIC POTENTIAL BLOCKS #
##################################
potential_blocks() {

  workdir=g_H_potential_blocks
  mkdir -p $workdir
  cd $workdir

  binsize=0.001
  blocksize=100000 # 100 ns
  ref_group="Lipids" # Reference group
  groups_water=("Water" "System") # groups which include water
  group_list=("CO" "CHO" "POPC" "DPPC" "Protein" "NA" "CL" "POPC_N" "POPC_P" "DPPC_N" "DPPC_P")

  # build list of groups which do not include water
  groups_nw=()
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      groups_nw+=("$g")
    fi
  done

  b=1
  while [[ $b -lt $lastframe ]]; do
    let e=${b}+${blocksize}-1
    mkdir ${b}-${e}
    cd ${b}-${e}

    # calculations using trajectories WITHOUT water
    for group in ${groups_nw[@]}
    do
      mkdir -p $group
      cd $group
      select="\"com of group $ref_group pbc; group $group\""
      echo $group

      # g_H_potential
      sem -j $maxjobs_nw g_H_potential -f $traj_nw -n $index_nw -s $structure_nw  -b $b -e $e -geo Radial -bin_size $binsize -select "$select" -dt $dt 
      
      cd ..
    done

    # calculations using trajectories WITH water
    for group in ${groups_water[@]}
    do
      mkdir -p $group
      cd $group
      select="\"com of group $ref_group pbc; group $group\""

      # g_H_potential
      sem -j $maxjobs g_H_potential -f $traj -n $index -s $structure  -b $b -e $e -geo Radial -bin_size $binsize -select "$select" -dt $dt 
      cd ..
    done

    cd ..
    let b=${b}+${blocksize}
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

  ref_group="Lipids" # Reference group
  group="Water" # water group

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
    echo "$ref_group $group" | sem -j $maxjobs g_sorient -f $traj -n $index -s $structure -b $begin -cbin $cbin -rbin $rbin -com -rmin $rmin -rmax $rmax -dt $dt 

    #next slice:
    rmin=$rmax
    rmax=$(echo "$rmin+$rstep" | bc -l)
    cd ..

  done

  cd ..

}


############################
# WATER ORIENTATION BLOCKS #
############################
sorient_blocks() {

  if ! [[  $lastframe ]]; then
    echo "-e not defined"
    return 2
  fi

  if ! [[  $radius ]]; then
    echo "-r not defined"
    return 2
  fi

  workdir=g_sorient_blocks
  mkdir -p $workdir
  cd $workdir

  ref_group="Lipids" # Reference group
  group="Water" # water group
  blocksize=100000 # 100 ns

  # Calculate sorient for 0.2nm slices
  rstep=0.1
  rmin=$(echo "$radius-$rstep" | bc -l)
  rmax=$(echo "$radius+$rstep" | bc -l)
  cbin=0.05
  rbin=0.05

  b=1
  while [[ $b -lt $lastframe ]]; do
    let e=${b}+${blocksize}-1
    mkdir -p ${b}-${e}
    cd ${b}-${e}

    # g_sorient
    echo "$ref_group $group" | sem -j $maxjobs g_sorient -f $traj -n $index -s $structure -b $b -e $e -cbin $cbin -rbin $rbin -com -rmin $rmin -rmax $rmax -dt $dt 

    cd ..
    let b=${b}+${blocksize}
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
  group_list=(HDL Protein Lipids Monolayer CO CHO POPC DPPC POPC_Protein DPPC_Protein)

  # Ion contacts
  ref_groups=()
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      ref_groups+=("$g")
    fi
  done
  if [[ -n $traj_nw ]]; then
    groups="NA CL"
    for ref_group in ${ref_groups[@]}; do
      # g_mindist
      echo "$ref_group $groups" | sem -j $maxjobs_nw g_mindist -f $traj_nw -n $index_nw -s $structure_nw -group -ng 2 -dt $dt -od "${ref_group}-ions_mindist.xvg" -on "${ref_group}-ions_numcount.xvg" -d $dist
    done
  fi

  # Water contacts
  ref_groups=()
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      ref_groups+=("$g")
    fi
  done
  if [[ -n $traj ]]; then
    groups="Water"
    for ref_group in ${ref_groups[@]}; do
      # g_mindist
      echo "$ref_group $groups" | sem -j $maxjobs g_mindist -f $traj -n $index -s $structure -group -ng 1 -dt $dt -od "${ref_group}-water_mindist.xvg" -on "${ref_group}-water_numcount.xvg" -d $dist 
    done
  fi


  cd ..
}



########
# SASA #
########
sas() {

  workdir=g_sas
  mkdir -p $workdir
  cd $workdir

  group_list=(CO CHO POPC DPPC Protein Monolayer Lipids HDL)
  ref_group="HDL"

  # Build target group list
  groups=()
  for g in ${group_list[@]}; do
    if [[ $(grep " $g " $index_nw) ]]; then
      groups+=("$g")
    fi
  done

  for group in  ${groups[@]}; do
    # g_sas
    echo "$ref_group $group" | sem -j $maxjobs_nw g_sas -f $traj_nw -n $index_nw -s $structure_nw -o $group-area.xvg -or $group-resarea.xvg -oa $group-atomarea.xvg -tv $group-volume.xvg -q $group-connelly.pdb -dt $dt 
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

  chains=(Chain_1 Chain_2 Chain_3 Chain_4 Protein)

  for chain in ${chains[@]}; do
    mkdir -p $chain
    cd $chain

    # do_dssp
    echo "$chain" | sem -j $maxjobs_nw do_dssp -f $traj_nw -n $index_nw -s $structure_nw -dt $dt 

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
    echo $group | sem -j $maxjobs_nw g_gyrate -f $traj_nw -n $index_nw -s $structure_nw -dt $dt -o ${group}.xvg
  done

  cd ..
}



#####################
# run main function #
#####################
main

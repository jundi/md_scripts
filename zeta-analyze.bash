#!/bin/bash

##############
# trajectory #
##############
traj="traj.xtc"
structure="topol.tpr"
begin=0



####################
# index.ndx groups #
####################
index="index.ndx"
declare -A group
group["System"]=0
group["CO"]=13
group["POPC"]=14
group["protein"]=1
group["NA"]=15
group["CL"]=16
group["Water"]=22
group["Lipids"]=26
group["HDL"]=27
group["POPC_P"]=28
group["POPC_N"]=29
group["POPC_Protein"]=30




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
      echo "rtfm"
      exit 1
  esac
  shift       
done



########
# main #
########
main() {
  mindist
}



#######
# RDF #
#######
rdf() {

  workdir=g_rdf
  mkdir -p $workdir
  cd $workdir

  ref_group="${group["lipids"]}"
  groups="${group["CO"]} ${group["POPC"]} ${group["Protein"]} ${group["Na"]} ${group["Cl"]} ${group["water"]} ${group["P"]} ${group["N"]}"

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

  # Reference group
  ref_group="${group["lipids"]}"
  groups="${group["CO"]} ${group["POPC"]} ${group["Protein"]} ${group["lipids"]} ${group["HDL"]} "

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
  ref_group="${group["lipids"]}"

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
  source ~/.local/gromacs-fftpack/bin/GMXRC.bash
  binsize=0.05
  # Reference group
  ref_group="${group["lipids"]}"
  group_names=( CO POPC protein Na Cl water system N P )

  for g in "${group_names[@]}"
  do
    mkdir $g
    cd $g
    calc_group="${group["$g"]}"
    select="com of group $ref_group pbc; group $calc_group"
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
  ref_group="${group["lipids"]}"
  # water group
  calc_group="${group["water"]}"

  # Calculate sorient for 0.5nm slices
  rstep=0.5
  rmin=0
  rmax=$rstep
  rmaxmax=10

  while [[ $(echo "$rmax < $rmaxmax" | bc -l) == 1 ]]; do # bash can't compare floats...

    # Create new directory for every slice 
    mkdir "$rmin-$rmax"
    cd "$rmin-$rmax"

    echo "$ref_group $calc_group" | g_sorient -f ../../$traj -n ../../$index -s ../../$structure -b $begin -com -rmin $rmin -rmax $rmax
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
  ref_group_names=( CO POPC protein lipids HDL POPC_protein )
  groups="${group["Na"]} ${group["Cl"]} ${group["water"]}"

  for rg in "${ref_group_names[@]}"; do
    ref_group="${group["$rg"]}"
    echo "$ref_group $groups" | g_mindist -f ../$traj -n ../$index -s ../$structure -b $begin -group -ng 3 -od "$rg-mindist.xvg" -on "$rg-numcount.xvg" -d $dist
  done

  cd ..
}



#####################
# run main function #
#####################
main

#!/bin/bash
#
# Script for sending multiple jobs to SLURM resource manager.
#

# The jobname of each job will be "$prefix$workdir".
prefix='7BP_' 
# Time limit for each job 
hours='72'
# The batch script to be used. 
batchscript='batch-parallel-taito.sh'   
# List of workdirs.
dirs=("n3.6028" "n3.6991" "n3.8008" "n3.8999") 



# get list of commands to be sent
cmds=()
for d in ${dirs[@]}; do

  workdir=$(readlink -f $d)
  lastjobid="0"
  dependency=""

  # check if there is job running in this directory.
  while read rid rdir; do

    rdir=$(readlink -f $rdir)
    if [[ $workdir == $rdir ]]; then
      echo "There is a job already running in $d."
      if [[ $rid -gt $lastjobid ]]; then
        dependency="-d afterok:${rid} "
        lastjobid=$rid
      fi
    fi

  done < <(squeue -u $USER -o "%i %Z")


  # build command
  cmd="sbatch -D ${d} -t ${hours}:00:00 -J ${prefix}${d} ${dependency}${batchscript} ${hours}"

  # print and save command
  echo "$cmd"
  cmds+=("$cmd")

done



# verify and send commands
echo ""
read -p "Send jobs? [yes|no]: " yn
case $yn in
  "yes")
  for cmd in "${cmds[@]}"; do
    eval "$cmd"
  done
  ;;
  *)
  echo "Canceling..."
  ;;
esac

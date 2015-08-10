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
# Check if there is already a job running, and start new job after it has
# finished.
deps='yes' # yes|no
# List of workdirs.
dirs=("n3.6028" "n3.6991" "n3.8008" "n3.8999") 



# get list of commands to be sent
cmds=()
for d in ${dirs[@]}; do

  # dependency flag
  if [[ $deps == "yes" ]];then
    jobid=$(squeue -u $USER -o "%i %j " | grep "${prefix}${d} " | cut -d " " -f 1)
    if [[ -n $jobid ]]; then
      dependency="-d afterok:${jobid} "
    else 
      dependency=""
    fi
  else
    dependency=""
  fi

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

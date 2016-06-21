#!/bin/bash
#
# Check the performance of Gromacs simulation from md.log files produced by
# mdrun.
#
# Usage: scaling.sh file1.log [file2.log]...
#

for f in "$@"; do
  cpus=$(grep "Using.*processes" $f | cut -d " " -f 2 | tail -n 1)
  nspday=$(grep "Performance" $f | tr -s " " | cut -d " " -f 2 | tail -n 1)
  if [[ -z "$nspday" ]]; then
    echo "$f: no stats found"
    continue
  fi
  nsphour=$(echo "scale=3; $nspday/24" | bc -l)
  nsphourpcpu=$(echo "scale=3; $nsphour/$cpus" | bc -l)
  nspdaypcpu=$(echo "scale=3; $nspday/$cpus" | bc -l)
  printf "%-12s %4s CPUS, %7s ns/day, %6s ns/hour, %6s ns/hour/cpu, %6s ns/day/cpu\n" "$f:" "$cpus" "$nspday" "$nsphour" "$nsphourpcpu" "$nspdaypcpu"
done

#!/bin/bash

for f in "$@"; do
  cpus=$(grep "Using.*processes" $f | cut -d " " -f 2 | tail -n 1)
  nspday=$(grep "Performance" $f | tr -s " " | cut -d " " -f 2 | tail -n 1)
  if [[ -z "$nspday" ]]; then
    continue
  fi
  nsphour=$(echo "scale=3; $nspday/24" | bc -l)
  nsphourpcpu=$(echo "scale=3; $nsphour/$cpus" | bc -l)
  nspdaypcpu=$(echo "scale=3; $nspday/$cpus" | bc -l)
  echo "$f $cpus CPUS,  $nspday ns/day, $nsphour ns/hour, $nsphourpcpu ns/hour/cpu, $nspdaypcpu ns/day/cpu"
done

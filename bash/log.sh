#!/bin/bash
# This script is just a shortcut to open md.log files produced GROMACS mdrun
# (Yes, I am that lazy). The advantage compared to using alias is that this
# also checks if there is multiple .log-files around...

ASK="Select file to view"
PS3="File: "
NOFILE="No log files found"
VIEWER="less +G"

maara=`ls *.log | wc -l`

if [[ $maara -gt 1 ]]; then
  echo $ASK
  select tiedosto in *log all; do
    case $tiedosto in
      *log) $VIEWER $tiedosto
      ;;
      all) $VIEWER *log
      ;;
      *) break
      ;;
    esac
    break
  done
elif [[ $maara -eq 1 ]]; then
  $VIEWER *log
else
  echo $NOFILE
fi

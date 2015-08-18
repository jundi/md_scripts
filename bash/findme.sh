#!/bin/bash
#
# Highlights lines which contain word "$USER"
#

# colors
highlight_color="\e[1;32m"
default_color="\e[0m"

# std or file?
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

# print and highlight
cat $input | while read words; do
  if [[ $words == *$USER* ]]; then
    echo -e "${highlight_color}${words}${default_color}"
  else
    echo -e "$words"
  fi
done

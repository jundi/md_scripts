#!/usr/bin/env python
#
# Rips one column from multicolumn xvg to new file. New file is named as
# "oldfile-column_n.xvg" if output file is not given.
#
# Example usage:
# $ python xvg_choose_column.py -i old_file.xvg -c 2 -o new_file.xvg
#

import argparse



# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True, nargs='+')
parser.add_argument('-o', required=False)
parser.add_argument('-c', required=True)
args = parser.parse_args()

infile_name = args.f
if args.o:
    outfile_suffix = args.o
else:
    outfile_suffix = '-column_'
column_num = args.c




# open files
for f in infile_name:
    infile = open(f,'r')
    outfile_name = '.'.join(f.split('.')[:-1]) + outfile_suffix + column_num + '.' + f.split('.')[-1]
    outfile = open(outfile_name,'w')

    for line in infile:
        if not (line.startswith('#') or line.startswith('@') or line.startswith('&')):
            words=line.split()
            print(words)
            outfile.write(words[0] + ' ' + words[int(column_num)] + '\n')

    # close files
    infile.close() 
    outfile.close()

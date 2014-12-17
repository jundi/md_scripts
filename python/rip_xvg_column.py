# Rips one column from multicolumn xvg to new file. New file is named as
# "oldfile-column_n.xvg"
#
# Example usage:
# $ python rip_xvg_column.py old_file.xv 2
#

import sys

infile_name = sys.argv[1]
column_num = sys.argv[2]
outfile_name = '.'.join(infile_name.split('.')[:-1]) + '-column_' + column_num + '.' + infile_name.split('.')[-1]


# open files
infile = open(infile_name,'r')
outfile = open(outfile_name,'w')

for line in infile:
	if line.startswith('#') or line.startswith('@'):
		outfile.write(line)
	else:
		words=line.split()
		outfile.write(words[0] + ' ' + words[int(column_num)] + '\n')

# close files
infile.close() 
outfile.close()

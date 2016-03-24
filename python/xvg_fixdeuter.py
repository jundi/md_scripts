#!/usr/bin/env python

'''Combines saturated and unsaturated deuter.xvg -plots, and fixes atom
numbering.

For example POPC oleyl tail:
xvg_fixdeuter -f saturated.xvg -u unsaturated.xvg -a 9 10 -o fidex.xvg
'''


import argparse

# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', required=True)
parser.add_argument('-u' )
parser.add_argument('-a', nargs='+')
parser.add_argument('-o', required=True)
args = parser.parse_args()

filename = args.f
unsat_filename = args.u
unsat_atoms = args.a
output_filename = args.o



# read saturated tail
tail=[]
with open(filename,'r') as file:
    for line in file:
        if line.startswith('#') or line.startswith('@'):
            continue
        else:
            tail.append(line.split()[1:])



# read unsaturated atoms
if unsat_filename:
    print("Tail is unsaturated. Replacing atoms {}.".format(", ".join(unsat_atoms)))
    with open(unsat_filename,'r') as file:
        itr_a=0
        for line in file:
            if line.startswith('#') or line.startswith('@'):
                continue
            else:
                atom=unsat_atoms[itr_a]
                tail[int(atom)-2] = line.split()[1:]
                itr_a = itr_a+1
else:
    print("Tail is saturated.")


# write file
itr=2
with open(output_filename,'w') as output:
    for t in tail:
        t_formatted = ''.join([x.rjust(15) for x in t])
        line = '{:>5} {}'.format(str(itr), t_formatted)
        output.write(line + "\n")
        itr=itr+1

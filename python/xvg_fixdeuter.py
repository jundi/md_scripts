#!/usr/bin/env python

'''Combines deuter.xvg -plots.'''

import argparse

# parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', required=True, nargs='+')
parser.add_argument('-u', nargs='+')
parser.add_argument('-a', nargs='+')
parser.add_argument('-o', required=True)
args = parser.parse_args()

filenames = args.s
unsat_filenames = args.u
unsat_atoms = args.a
output_filename = args.o



# read saturated tails
tails=[]
for filename in filenames:
    with open(filename,'r') as file:
        tail=[]
        for line in file:
            if line.startswith('#') or line.startswith('@'):
                continue
            else:
               tail.append(line.split()[1])
        tails.append(tail)



# read unsaturated atoms
itr=0
for filename in unsat_filenames:
    if filename:
        atoms=unsat_atoms[itr].split()
        print("sn{}-tail is unsaturated. Replacing atoms {}.".format(itr+1, ", ".join(atoms)))
        with open(filename,'r') as file:
            itr_a=0
            tail=[]
            for line in file:
                if line.startswith('#') or line.startswith('@'):
                    continue
                else:
                    atom=atoms[itr_a]
                    tails[itr][int(atom)-2] = line.split()[1]
                    itr_a = itr_a+1
    else:
        print("sn{}-tail is saturated.".format(itr+1))

    itr=itr+1



# length of longest tail
maxlength=0
for t in tails:
    maxlength = max(maxlength, len(t))

# write file
itr=0
with open(output_filename,'w') as output:
    while itr < maxlength:
        words=[str(itr+2)]
        for t in tails:
            if itr < len(t):
                words.append(t[itr])
            else:
                words.append("")
        words=[w.rjust(15) for w in words]
        line = "".join(words)
        output.write(line + "\n")
        itr=itr+1

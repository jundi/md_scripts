#!/bin/python
'''
This script names data sets of .agr-file according to comments.
'''
#import sys
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-c', default=-1, help='which word to use')
parser.add_argument('-f', help='input')
parser.add_argument('-o', help='output')
args = parser.parse_args()

if args.f:
    infile_name=args.f
else:
    print('input file missing')
    exit()
if args.o:
    outfile_name=args.o
else:
    print('output file missing')
    exit()



def choose_legend(comment,n):
    words=re.split('\.|/|"',comment)
    nonempty_words=[x for x in words if x]
    if n == -1:
        print(nonempty_words)
        n=int(input('Which? (0-' + str(len(nonempty_words)-1) + ')'))
    return (nonempty_words[n],n)



newfile=[]
comment_num=int(args.c)
with open(infile_name,'r') as infile:
    for line in infile:
        words=line.split()
        if len(words) > 2 and words[2] == 'comment':
            comment=words[3]
        if len(words) > 2 and words[2] == 'legend':
            legend, comment_num = choose_legend(comment,comment_num)
            print(legend)
            newfile.append(' '.join(words[0:3]) + ' "' + legend+ '"' + '\n')
        else:
            newfile.append(line)


with open(outfile_name,'w') as outfile:
    for line in newfile:
        outfile.write(line)


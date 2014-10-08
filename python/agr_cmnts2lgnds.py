#!/bin/python
'''
This script names data sets of .agr-file according to comments.
'''
import sys
import re

infile_name=sys.argv[1]
outfile_name=sys.argv[2]

def choose_legend(comment,n):
    words=re.split('\.|/|"',comment)
    nonempty_words=[x for x in words if x]
    if n == -1:
        print(nonempty_words)
        n=int(input('Which? (0-' + str(len(nonempty_words)) + ')'))
    return (nonempty_words[n],n)



newfile=[]
comment_num=-1
with open(infile_name,'r') as infile:
    for line in infile:
        words=line.split()
        if len(words) > 2 and words[2] == 'comment':
            comment=words[3]
        if len(words) > 2 and words[2] == 'legend':
            legend, comment_num = choose_legend(comment,comment_num)
            print(legend)
            newfile.append(' '.join(words[0:3]) + legend + '\n')
        else:
            newfile.append(line)


with open(outfile_name,'w') as outfile:
    for line in newfile:
        outfile.write(line)


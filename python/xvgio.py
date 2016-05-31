#!/bin/python

import numpy

def read(filename):
    comment = []
    data = []
    with open(filename,'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith('@'):
                comment.append(line)
            else:
                points = [float(x) for x in line.split()]
                data.append(points)
    return (numpy.array(data),comment)


def write(filename,data,comment):
    with open(filename,'w') as file:
        for row in comment:
            file.write(row)
        for row in data:
            line=''
            for word in row:
                line=line + str('%11e' % word) + ' '
            file.write(line + '\n')

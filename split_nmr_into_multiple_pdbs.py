#! /usr/bin/python2
# Author Firas Khatib
from os import popen
import sys

input = open(sys.argv[1], 'r').readlines()

model = -1
for i in input:
    line = i.split()
    if len(line) >= 1:
        if line[0] == "MODEL":
            model = model + 1 # int(line[1])
            output = open(sys.argv[1]+".model."+str(model)+".pdb", 'w')

        if model >= 0:
            output.write(i)

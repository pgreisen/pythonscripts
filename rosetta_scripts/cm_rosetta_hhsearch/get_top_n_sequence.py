#!/usr/bin/env python
import sys,os
import argparse

hhsearchfile = sys.argv[1]
no_of_alignments = int(sys.argv[2])

newfile = []

with open(hhsearchfile, 'r') as f:
    for line in f:
        tmpline = line.split()
        print(tmpline)
        if( len(tmpline) == 2 and  tmpline[0] == 'No' and int(tmpline[1]) > no_of_alignments ):
            break
        else:
            newfile.append( line )

with open("reduced_hhsearch.hhr",'w') as f:
    for line in newfile:
        f.write(line)


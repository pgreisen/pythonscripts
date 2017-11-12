from os import popen
import sys

input = open(sys.argv[1], 'r').readlines()

with open("new.pdb",'w') as f:
    for line in input:
        # tmp = line.split()
        if len(line) > 4 and line[0:4] == "HETA":

            if line[13:14] == "H":
                continue
            elif line[13:16] == "C21":
                continue
            elif line[13:16] == "C22":
                continue
            else:
                f.write(line)
        else:
            f.write(line)


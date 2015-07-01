from os import popen
import sys

histag = ["H28"]

input = open(sys.argv[1], 'r').readlines()

with open("new_file.pdb",'w') as f:
    for line in input:
        tmp = line.split()
        if len(line) > 4 and line[0:4] == "ATOM" or line[0:4] == "HETA":

            if tmp[2] == "H28":
                continue
            else:
                f.write(line)
        else:
            f.write(line)


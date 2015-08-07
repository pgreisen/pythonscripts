from os import popen
import sys

input = open(sys.argv[1], 'r').readlines()

with open("new.pdb",'w') as f:
    for line in input:
        # tmp = line.split()
        if len(line) > 4 and line[0:4] == "HETA":

            if line[13:14] == "H":
                continue
            elif line[13:16] == "C6 ":
                newline = line[0:13]+str("C10")+line[16:]
                f.write(newline)
            elif line[13:16] == "C7 ":
                newline = line[0:13]+str("C6 ")+line[16:]
                f.write(newline)

            elif line[13:16] == "C8 ":
                newline = line[0:13]+str("C7 ")+line[16:]
                f.write(newline)

            else:
                f.write(line)
        else:
            f.write(line)


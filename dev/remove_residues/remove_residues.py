from os import popen
import sys

histag = ["68","69","70","71","72","73","74","75","76","77"]

input = open(sys.argv[1], 'r').readlines()

model = -1
for line in input:
    tmp = line.split()
    if len(tmp) >= 1:
        if( tmp[0] == "MODEL" ): #or tmp[0] == "TER":
            model = model + 1 # int(line[1])
            output = open(sys.argv[1]+".model."+str(model)+".pdb", 'w')
        if( model >= 0 and len(line) > 59 and line[0:4] == "ATOM" ):
            if( str(line[22:26]).strip() not in histag):
                output.write(line)

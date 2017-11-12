import os,sys

f1 = sys.argv[1]

positions = []

with open(f1,'r') as f:
    for line in f:
        tmp = line.split()
        for i in tmp:
            positions.append( i )

with open("position_patchdock",'w') as pd:
    for character in positions:
        pd.write( character +" A \n")

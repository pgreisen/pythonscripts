import sys
'''
Takes a text file with mutations and create a single object for each of them in pymol.

'''
filename = sys.argv[1]
method = sys.argv[2]

pymollist = []
prefix = "create "
postfix = "_"+method+" , resi "

with open(filename,'r') as f:
    for line in f:
        tmpobj = prefix+line.strip()+postfix+str(line.strip()[1:-1])
        pymollist.append(tmpobj)

with open("mutations.pml",'w') as f:
    for line in pymollist:
        f.write(line+"\n")


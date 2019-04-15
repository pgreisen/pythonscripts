import sys

def get_dictionary(filename):
    sequences = {}
    with open(filename,'r') as f:
        for line in f:
            if(line[0] == ">"):
                tmp = line.split()
                ids = tmp[0]
                sequences[ids] = ""
            else:
                tmp = line.strip()
                sequences[ids] = sequences[ids] + tmp
    return sequences

dfile1 = sys.argv[1]
designs = get_dictionary(dfile1)

dfile2 = sys.argv[2]
alldesigns = get_dictionary(dfile2)

newround = {}
for key in alldesigns.keys():
    if(key not in designs.keys()):
        newround[key] = alldesigns[key]
    else:
        print(key)
        
with open("nextround.fasta",'w') as f:
    for key in newround.keys():
        f.write(key+"\n")
        f.write(newround[key]+"\n")

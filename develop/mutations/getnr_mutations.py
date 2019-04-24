import sys
import numpy as np
import matplotlib.pyplot as plt
nrmuts = []
inputfile = sys.argv[1]

nrmt = {}


with open(inputfile, 'r') as f:
    for line in f:
        if(line[0] == ">"):
            tmp = len(line.split('_'))
            nrmuts.append(tmp)
            tmp_key = "Muts_"+str(tmp)
            if(tmp_key not in nrmt.keys()):
                nrmt[tmp_key] = 0
            nrmt[tmp_key] += 1
                

plt.hist(nrmuts)
print(np.mean(nrmuts) )
plt.show()

for key in nrmt.keys():
    print(key, str(nrmt[key]) )

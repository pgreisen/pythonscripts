import sys
import numpy as np
import matplotlib.pyplot as plt
nrmuts = []
inputfile = sys.argv[1]

with open(inputfile, 'r') as f:
    for line in f:
        if(line[0] == ">"):
            tmp = len(line.split('_'))-1
            nrmuts.append(tmp)

plt.hist(nrmuts)
print(np.mean(nrmuts) )
plt.show()
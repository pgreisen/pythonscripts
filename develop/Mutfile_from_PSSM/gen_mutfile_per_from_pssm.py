import sys

NUMBER_OF_THREADS = 5


'''
# generate mutation file for ddG cartesian
header="total 1\n1\n"
for i,j,k in zip(dftot.AA_WT, dftot.AA_Pos,dftot.AA):
    with open(i+str(j)+k+".mutfile",'w') as f:
        f.write(header)
        f.write(i+" "+str(j)+" "+k)
'''


def get_mutfile_pssm(mutfile):
    seq = ""
    with open(fastafile,'r') as f:
        for line in f:
            if(line[0] != ">" ):
                seq += line.strip()
    return seq


mutfile = sys.argv[1]


# dictionary to hold data
muts = {}
dummy_key = 0

with open(mutfile,'r') as f:
    for line in f:
        if(dummy_key % NUMBER_OF_THREADS == 0):
            tmpkey = str(dummy_key)
            muts[tmpkey] = []
        else:
            wt_ = line[0]
            aa_ = line.strip()[-1]
            pos_ = line.strip()[1:-1]
            tmp_ = "1\n"+wt_+" "+str(pos_)+" "+aa_+"\n"
            muts[tmpkey].append(tmp_)
        dummy_key += 1

# write to file
for key in muts.keys():
    with open("muts_"+str(key)+".mutfile",'w') as f:
        lnght = len(muts[key])
        f.write("total "+str(lnght)+"\n")
        for i in muts[key]:
            f.write(i)
from itertools import *

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def get_cstfile(cstfilename,lst):
    with open(cstfile,'r') as f:
        for line in f:
            lst.append(line)
    return lst


def write_cstfile(newnamecst,newcstfile):
    with open("cstfiles/"+newnamecst+".cst",'w') as f: 
        for line in newcstfile:
            f.write(line)


a = ["PP_1.cst", "PP_2.cst", "PP_5.cst", "TRP_3.cst", "tstack_1.cst", "tstack_4.cst","PP_3.cst","TRP_1.cst","TRP_4.cst","tstack_2.cst","tstack_5.cst","PP_4.cst","TRP_2.cst","TRP_5.cst","tstack_3.cst"]

b = combinations(a,3)

cstfiles = []

for j in b:
    lst= list(j)

    key = []
    for i in lst:
        tmp = i.split('_')
        key.append(int (tmp[1][0]) )

    if(len( set(key) ) == 3):
        cstfiles.append(lst)


for block in cstfiles:
    newcstfile = []
    newnamecst = ""
    for cstfile in block:
        newnamecst = newnamecst + cstfile.split(".")[0]+"_"
        newcstfile = get_cstfile(cstfile, newcstfile)
    write_cstfile(newnamecst,newcstfile)
        
        


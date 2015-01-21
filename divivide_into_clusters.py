import sys,os,subprocess,shutil

cluster_file = sys.argv[1]

with open(cluster_file,'r') as f:
    for line in f:
        if line[0] == '>':
            tmp = line.split()
            cluster = tmp[0]+tmp[1]
            print "The cluster is: ", cluster[1:]
            os.mkdir(cluster[1:])
            
        else:
            tmp = line.split()
            pdbid = ""
            for j in str(tmp[2]):
                if(j == '>'):
                    continue
                elif( j =='.'):
                    continue
                elif( j =='A'):
                    continue
                else:
                    pdbid = pdbid + j
            pdbid = pdbid+".pdb"
            print pdbid+" into the following directory: "+cluster[1:]

            shutil.copy(pdbid, cluster[1:])


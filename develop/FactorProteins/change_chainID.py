import os,shutil,sys,subprocess

'''
change the chain id for protein chain

'''

def main():
    
    pdbfile = sys.argv[1]
    
    newchainid = "A" #sys.argv[2]

    filename = open(pdbfile,'r')

    tmpfile = open('tmp.pdb','w')
    
    for line in filename:

        #        if line[0:4] == 'HETA':
        if line[0:4] == "ATOM":
            newline = str(line[0:20])+' '+str(newchainid)+str(line[22:])
            newline = newline[0:26]+' '+newline[27:]
            tmpfile.write(newline)
        else:
            tmpfile.write(line)
    tmpfile.close()
    filename.close()

    mvfl = 'mv tmp.pdb '+str(pdbfile)

    subprocess.Popen(mvfl,shell=True).wait()

    
if __name__ == "__main__":
    main()


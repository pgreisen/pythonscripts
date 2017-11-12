import os,shutil,sys,subprocess

'''
change the chain id for protein chain

'''

def main():
    
    pdbfile = sys.argv[1]
    
    new_resname = sys.argv[2]

    filename = open(pdbfile,'r')

    tmpfile = open('tmp.pdb','w')
    
    for line in filename:
        newline = str(line[0:17])+new_resname+str(line[20:])
        tmpfile.write(newline)

    tmpfile.close()
    filename.close()

    mvfl = 'mv tmp.pdb '+str(pdbfile)

    subprocess.Popen(mvfl,shell=True).wait()

    
if __name__ == "__main__":
    main()


import os,shutil,sys,subprocess

def main():
    
    pdbfile = sys.argv[1]
    
    filename = open(pdbfile,'r')

    tmpfile = open('tmptmptmp.pdb','w')
    
    for line in filename:

        if line[0:4] == 'ATOM':
            newline = str(line[0:20])+' A'+str(line[22:])
            tmpfile.write(newline)
        else:
            tmpfile.write(line)
    tmpfile.close()
    filename.close()

    mvfl = 'mv tmptmptmp.pdb '+str(pdbfile)

    subprocess.Popen(mvfl,shell=True).wait()

    
if __name__ == "__main__":
    main()


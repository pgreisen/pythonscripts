#!/usr/bin/env python2.7
import sys

'''
Take a pdbfile and converts it into a alanine scaffold.

'''


def readfile(pdbfile):
    tmp_file = open(pdbfile,'r')
    pfile = tmp_file.readlines()
    tmp_file.close()
    return pfile


def insert_alanine(pfile,pname):
    pwrite = open(pname+"_alanine_scaffold.pdb",'w')
    for line in pfile:
        if(line[0:4] == "ATOM"):
            print line[12:15].strip()
            if(line[12:15].strip() in ['N','CA','C','O','CB']):
                new_line = line[0:17]+"ALA"+line[20:]
                pwrite.write(new_line)
    pwrite.close()
    


def main():
    # pdbfile
    pdbfile = sys.argv[1]
    pdbname = pdbfile[0:4]
    pfile = readfile(pdbfile)
    
    insert_alanine(pfile,pdbname)

    
if __name__ == "__main__":
    main()

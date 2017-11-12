#! /usr/bin/python2
from os import popen
import sys, shutil, os, subprocess, argparse
from collections import defaultdict


class SplitMatches:

    def __init__(self):
        self.store_pdbs = {}
        self.matches = defaultdict(list)
        self.models = []
        self.pdbfile = ""
        self.template_model = "MODEL    1"


    def replace_residues(self):
        '''
        here we replace the different residues
        '''


    def get_models(self):
        '''


        '''
        for key in self.matches:
            if key == self.template_model:
                with open(key+".pdb",'w') as f:
                    for line in self.matches[key]:
                        f.write(line)
            else:
                with open(key+".pdb",'w') as f:

                    for ln1 in self.matches[self.template_model]:
                        boolline = 0
                        for ln2 in self.matches[key]:
                            tmp1 = ln1.split()
                            tmp2 = ln2.split()
                            if( len(tmp1) > 5 and len(tmp2) > 5):

                                atm = tmp1[2]
                                resn = tmp1[3]
                                resi = tmp1[5]

                                tmp2 = ln2.split()
                                atm2 = tmp2[2]
                                resn2 = tmp2[3]
                                resi2 = tmp2[4]

                                if(atm == atm2 and resn == resn2 and resi == resi2):
                                    print atm,atm2,resn,resn2,resi,resi2

                                if(atm == atm2 and resn == resn2 and resi == resi2):
                                    # print ln1,ln2
                                    boolline =1
                                    newline = ln2
                        if( boolline == 0  ):
                            f.write(ln1)
                        else:
                            f.write( newline )




    def set_pdbs(self,pdbfile):
        '''
        '''
        key = ""
        for line in pdbfile:
            if(len(line) > 10 and line[0:5] == "MODEL"):
                key = line.strip()
                # debug
                print key
            if(len(key) > 1 and line[0:6] != "ENDML"):
                self.matches[key].append(line)


    def main(self):
        parser = argparse.ArgumentParser(description="Split matches into multiple pdbs")
        parser.add_argument('-f',dest="pdbfile", help='PDB file' )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        input = open(self.pdbfile, 'r').readlines()
        self.set_pdbs(input)

        self.get_models()




if __name__ == '__main__':
    run = SplitMatches()
    run.main()
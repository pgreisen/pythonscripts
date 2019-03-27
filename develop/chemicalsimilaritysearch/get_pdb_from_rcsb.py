#import urllib2
import urllib 
import gzip
import os,sys, argparse, subprocess
from CleanPDB import *
from collections import OrderedDict
'''

This script will search the PDB for a chemical fragment and return a match based on the similarity between them. It assumes
that the user has access to openbabel and especially the two executables: babel and obfit.


'''

class GetPDBs:

    def __init__(self):
        self.pdbname = ""

    def convert_to_pdb_gz(self,pdbid):
        '''                                                                                                                       
        # Requires PDB ID                                                                                                         
        # Returns gz-file name                                                                                                    
        @param pdbid:                                                                                                             
        @return:filename                                                                                                          
        '''
        pdbid = pdbid.lower()
        filename = "pdb"+pdbid+".ent.gz"

        return filename



    def writes_to_pdb_format(self,filename,pdbname):
        # Uncompress the archive, delete when done
        f = gzip.open(filename, 'rb')
        f_content = f.read()
        f.close()
        outfile = open(pdbname+'.pdb', 'wb')
        outfile.writelines(f_content)
        outfile.close()
        os.remove(filename)

    def get_gz_pdbfile(self,pdbname):
        '''
        # Requires Protein Data Bank url
        # Writes pdb file
        @param pdbname:
        @return:
        '''
        ##import Bio
        ##from Bio.PDB import PDBList
        ##'''Selecting structures from PDB'''
        ##pdbl = PDBList()
        #PDBlist2=['4B97','4IPH','4HNO','4HG7','4IRG','4G4W','4JKW','4IPC','2YPM','4KEI']
        #for i in PDBlist2:
        ##pdbl.retrieve_pdb_file(pdbname,pdir='PDB')
        ##'''
        gz_name = self.convert_to_pdb_gz(pdbname)
        print(gz_name)
        url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'+gz_name[4:6]+'/'+gz_name
        print(url)
        try:
            urllib.urlretrieve(url, gz_name)
            #self.writes_to_pdb_format(gz_name,pdbname)
        except:
            print("Exit")
        ###'''

    def main(self):

        parser = argparse.ArgumentParser(description="Takes a pdb with a chemical fragment and seraches the PDB for the occurance of this fragment.")
        # get the initial rosetta design as input
        parser.add_argument("-s", dest="pdbname", help="This file contains the coordinates of the chemical fragment" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.get_gz_pdbfile( self.pdbname )

if __name__ == "__main__":
    run = GetPDBs()
    run.main()

import urllib2
import urllib 
import gzip
import os,sys, argparse, subprocess
from collections import OrderedDict

'''

This script will download a pdb from the Protein Databank(PDB)

'''

class GetPDBfDatabase:

    def __init__(self):
        self.maxlength = "200"
        self.resolution = "2.0"
        self.similarity = "1.0"
        self.babelbin = "/work/greisen/ExternalProgram/openbabel-2.3.2/bin/"
        self.babel = "babel"
        self.obfit = "/work/greisen/ExternalProgram/bin/obfit"
        self.format = "pdb"
        self.results = []
        self.smi_directly = 0
        self.smiles = ""
        self.pdbname = ""


    def writes_to_pdb_format(self,filename,pdbname):
        # Uncompress the archive, delete when done
        f = gzip.open(filename, 'rb')
        f_content = f.read()
        f.close()
        outfile = open(pdbname+'.pdb', 'wb')
        outfile.writelines(f_content)
        outfile.close()
        os.remove(filename)


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




    def get_gz_pdbfile(self,pdbname):
        '''
        # Requires Protein Data Bank url
        # Writes pdb file
        @param pdbname:
        @return:
        '''
        gz_name = self.convert_to_pdb_gz(pdbname)
        url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'+gz_name[4:6]+'/'+gz_name
        try:
            urllib.urlretrieve(url, gz_name)    
            self.writes_to_pdb_format(gz_name,pdbname)

        except:
            print "Not possible to retrieve the following address: ",url



    def main(self):
        parser = argparse.ArgumentParser(description="Download pdb from PDB database.")
        parser.add_argument("-s", dest="pdbname", help="This file contains the coordinates of the chemical fragment" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])
        
        self.get_gz_pdbfile(self.pdbname)

if __name__ == "__main__":
    run = GetPDBfDatabase()
    run.main()

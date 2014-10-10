#!/usr/bin/python
from optparse import OptionParser
import urllib2
import urllib
import gzip
import os
import shutil
# from shutil import copy


class GetPDBTMalign:

    def __init__(self):
        self.pdbids = []
        self.threshold = 0.5

        self.chains = {}



    def get_pdbs(self,datafile):
        with open(datafile,'r') as f:
            for line in f:
                tmp = line.split()

                if( float(tmp[7]) > self.threshold ):
                    self.pdbids.append( str(tmp[1].split('/')[6])[0:4]  )


    def writes_to_pdb_format(self,filename,pdbname):
        # Uncompress the archive, delete when done
        f = gzip.open(filename, 'rb')
        f_content = f.read()
        f.close()

        outfile = open(pdbname+'.pdb', 'wb')
        outfile.writelines(f_content)
        outfile.close()

        os.remove(filename)

    # Requires PDB ID
    # Returns gz-file name
    def convert_to_pdb_gz(self,pdbid):
        pdbid = pdbid.lower()
        filename = "pdb"+pdbid+".ent.gz"
        return filename



    def get_chains_from_pdb(self,pdbfile):

        with open(pdbfile+".pdb",'r') as f:
            for line in f:

                if( line[0:6] == "COMPND"):
                    tmp = line.split()
                    if( tmp[2] == "CHAIN:"):
                        # print "split of line",tmp[3:]  #, line
                        self.chains[pdbfile] = tmp[3:]


    # Requires Protein Data Bank url
    # Writes pdb file
    def get_gz_pdbfile(self, pdbname):
        gz_name = self.convert_to_pdb_gz(pdbname)
        url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'+gz_name[4:6]+'/'+gz_name
        urllib.urlretrieve(url, gz_name)
        self.writes_to_pdb_format(gz_name,pdbname)


    def main(self):
        parser = OptionParser()
        parser.add_option('-f',dest='datafile',help='data file from TMalign')

        parser.add_option('-t',dest='threshold', default='0.5', help='Threshold for the TMalign default=0.5',type=float)

        (options, args) = parser.parse_args()

        assert isinstance(options.threshold, float)
        self.threshold = options.threshold

        self.get_pdbs(options.datafile)

        for pdbid in self.pdbids:

            self.get_gz_pdbfile( pdbid )
            self.get_chains_from_pdb( pdbid )


        if( 1 == 1):
            try:
                os.mkdir("Monomer")
            except:
                print "Monomers exists"
            try:
                os.mkdir("Dimer")
            except:
                print "Dimer exists"
            try:
                os.mkdir("Oligomers")
            except:
                print "Oligomers exists"

        with open("pdbfile.txt",'w') as f:

            for key in self.chains:

                value = len( self.chains[key] )

                pdbfile = key+".pdb"

                if( value == 1):
                    shutil.copy( pdbfile, "Monomer" )
                elif ( value == 2 ):
                    shutil.copy( pdbfile,"Dimer")
                else:
                    shutil.copy( pdbfile, "Oligomers")

                f.write(str(key)+" "+str(self.chains[key])+"\n")








if __name__ == "__main__":
    run = GetPDBTMalign()
    run.main()

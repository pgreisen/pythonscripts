import urllib2
import urllib 
import gzip
import os,sys, argparse, subprocess
from collections import OrderedDict
'''

This script will search the PDB for a chemical fragment and return a match based on the similarity between them. It assumes
that the user has access to openbabel and especially the two executables: babel and obfit.


'''

class SearchPDBwSmiles:


    def get_pdbs_from_PDB(self):
        '''
        downloads the pdbs from the PDB based on SMILES on ligands
        @return:
        '''

        url = 'http://www.rcsb.org/pdb/rest/search'

        queryText = """
<?xml version="1.0" encoding="UTF-8"?>

<orgPdbQuery>

<queryType>org.pdb.query.simple.SequenceQuery</queryType>

<description>Sequence Search (Structure:Chain = 4HHB:A, Expectation Value = 10.0, Search Tool = BLAST)</description>

<sequence>MTDRYSFSLTTFSPSGKLGQIDYALTAVKQGVTSLGIKATNGVVIATEKKSSSPLAMSETLSKVSLLTPDIGAVYSGMGP
DYRVLVDKSRKVAHTSYKRIYGEYPPTKLLVSEVAKIMQEATQSGGVRPFGVSLLIAGHDEFNGFSLYQVDPSGSYFPWK
ATAIGKGSVAAKTFLEKRWNDELELEDAIHIALLTLKESVEGEFNGDTIELAIIGDENPDLLGYTGIPTDKGPRFRKLTS
QEINDRLEAL</sequence>

<eCutOff>30.0</eCutOff>

<searchTool>blast</searchTool>

<sequenceIdentityCutoff>95</sequenceIdentityCutoff>

</orgPdbQuery>

"""
        print "querying PDB...\n"
        req = urllib2.Request(url, data=queryText)

        f = urllib2.urlopen(req)

        # Contains pdb ids
        result = f.read().rstrip().split('\n')
        if len(result) > 1:

            print "Found number of PDB entries:", result.count('\n')

        else:

            print "Failed to retrieve results"


        for i in result:
            print "Downloading PDB structure '%s'..." % i
            #self.get_gz_pdbfile(i)

    def main(self):

        parser = argparse.ArgumentParser(description="Takes a pdb with a chemical fragment and seraches the PDB for the occurance of this fragment.")
        # get the initial rosetta design as input
        parser.add_argument("-s", dest="fragment_file", help="This file contains the coordinates of the chemical fragment" )
        parser.add_argument("--maxlength", dest="maxlength", help="The max length of the protein (Default=200 aa)", default="200", type=str )
        parser.add_argument("--resolution", dest="resolution", help="Resolution of crystal structure (Default=2.0)", type=str, default="2.0" )
        parser.add_argument("--similarity", dest="similarity", help="The chemical similarity between the fragment searched (Default=1.0)", type=str )
        #parser.add_argument("--babel", dest="babel", help="The path to the executable for openbabel (Default - the dig system in the Bakerlab" )
        parser.add_argument("--format", dest="format", help="Format to convert from (Default=pdb)",default="pdb" )
        parser.add_argument("--smi_directly", dest="smi_directly",default=0, help="The input file has to be query.smi to pass SMILE directly" )
        parser.add_argument("--smiles", dest="smiles",default=0, help="The input file has to be query.smi to pass SMILE directly" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])
        
        print self.similarity


        self.get_pdbs_from_PDB()


if __name__ == "__main__":
    run = SearchPDBwSmiles()
    run.main()

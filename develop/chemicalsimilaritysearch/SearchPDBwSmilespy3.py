import urllib 
import gzip
import os,sys, argparse, subprocess
from CleanPDB import *
from collections import OrderedDict
import requests
'''

This script will search the PDB for a chemical fragment and return a match based on the similarity between them. It assumes
that the user has access to openbabel and especially the two executables: babel and obfit.


'''

class SearchPDBwSmiles:

    def __init__(self):
        self.smiles = ""
        self.maxlength = "500"
        self.resolution = "2.5"
        self.similarity = "1.0"
        self.babelbin = "/work/greisen/ExternalProgram/openbabel-2.3.2/bin/"
        self.babel = "babel"
        self.obfit = "/work/greisen/ExternalProgram/bin/obfit"
        self.format = "pdb"
        self.results = []
        self.smi_directly = 0
        self.smiles = ""


    def get_query_pdb(self):
        '''
        downloads the pdbs from the PDB based on SMILES on ligands
        @return:
        '''

        queryText = """
<?xml version="1.0" encoding="UTF-8"?>

<orgPdbCompositeQuery version="1.0">

<queryRefinement>
<queryRefinementLevel>0</queryRefinementLevel> 

<orgPdbQuery>

<queryType>org.pdb.query.simple.ResolutionQuery</queryType>

<description>ResolutionQuery: refine.ls_d_res_high.comparator=between refine.ls_d_res_high.min=0.0 refine.ls_d_res_high.max="""+str(self.resolution)+""" </description>

<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>

<refine.ls_d_res_high.min>0.0</refine.ls_d_res_high.min>

<refine.ls_d_res_high.max>"""+str(self.resolution)+"""</refine.ls_d_res_high.max>

</orgPdbQuery> 

</queryRefinement>
  
<queryRefinement>                                                                                               
 
<queryRefinementLevel>1</queryRefinementLevel>  
<conjunctionType>and</conjunctionType> 

  <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
  <description>Chemical ID(s):  CL and Polymeric type is Any</description>
    <chemCompId>HEM</chemCompId>
    <polymericType>Any</polymericType>

</queryRefinement>

</orgPdbCompositeQuery>

"""
        return queryText



    def get_xtal_by_smile_search(self):
        
        url = 'http://www.rcsb.org/pdb/rest/search'

        query_text = self.get_query_pdb()
        
        print("Query: %s" % query_text)

        print("Querying RCSB PDB REST API...")

        header = {'Content-Type': 'application/x-www-form-urlencoded'}

        response = requests.post(url, data=query_text, headers=header)
        print(response, response.text)
        print("Debug")
        return response


        

    def set_smiles(self):
        with open("query.smi") as f:
            for line in f:
                self.smiles = line.split()[0]
                self.smiles.strip()

        print( "The PDB will be searched with this chemical fragment: "+self.smiles+"\n")

    def main(self):

        parser = argparse.ArgumentParser(description="Takes a pdb with a chemical fragment and seraches the PDB for the occurance of this fragment.")
        # get the initial rosetta design as input
        parser.add_argument("-s", dest="fragment_file", help="This file contains the coordinates of the chemical fragment" )
        parser.add_argument("--maxlength", dest="maxlength", help="The max length of the protein (Default=500 aa)", default="500", type=str )
        parser.add_argument("--resolution", dest="resolution", help="Resolution of crystal structure (Default=2.0)", type=str, default="2.0" )
        parser.add_argument("--similarity", dest="similarity", help="The chemical similarity between the fragment searched (Default=1.0)", type=str )
        #parser.add_argument("--babel", dest="babel", help="The path to the executable for openbabel (Default - the dig system in the Bakerlab" )
        parser.add_argument("--format", dest="format", help="Format to convert from (Default=pdb)",default="pdb" )
        parser.add_argument("--smi_directly", dest="smi_directly",default=0, help="The input file has to be query.smi to pass SMILE directly" )
        parser.add_argument("--smiles", dest="smiles",default=0, help="The input file has to be query.smi to pass SMILE directly" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])
        

        ##import pdb;pdb.set_trace()
        self.set_smiles()


        ##self.get_pdbs_from_PDB()
        xtal_search = self.get_xtal_by_smile_search()
        print(xtal_search.text)

        for pdbid in xtal_search.text.split("\n"):
            print(pdbid)
        #self.align_pdb_obfit()


if __name__ == "__main__":
    run = SearchPDBwSmiles()
    run.main()

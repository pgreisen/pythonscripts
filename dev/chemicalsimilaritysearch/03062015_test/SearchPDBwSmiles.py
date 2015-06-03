import urllib2
import urllib 
import gzip
import os,sys, argparse, subprocess
from CleanPDB import *
from collections import OrderedDict
'''

This script will search the PDB for a chemical fragment and return a match based on the similarity between them. It assumes
that the user has access to openbabel and especially the two executables: babel and obfit.


'''

class SearchPDBwSmiles:

    def __init__(self):
        # methoxy
        ##smiles='COC1=CC=CC=C1'
        self.smiles = ""
        self.maxlength = "200"
        self.maxresolution = "2.0"
        self.similarity = "1.0"
        self.babelbin = "/work/greisen/ExternalProgram/openbabel-2.3.2/bin/"
        self.babel = "babel"
        self.obfit = "/work/greisen/ExternalProgram/bin/obfit"
        self.format = "pdb"
        self.results = []


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

        urllib.urlretrieve(url, gz_name)
    
        self.writes_to_pdb_format(gz_name,pdbname)

        # clean the pdb
        # 22-05-2015
        cleanpdb = CleanPDB()

        pdbfile = pdbname+".pdb"

        pdb_by_chain, hetatoms = cleanpdb.get_chains(pdbfile)

        dummy = 1
        for key in pdb_by_chain:
            for key2 in  hetatoms:

                key2_split = key2.split('_')

                if(key == key2_split[1] ):
                    ## import pdb; pdb.set_trace()
                    # print "2########################################2",key,key2[0], key2
                    # with open(pdbname+"_chain"+key+".pdb",'w') as f:
                    ##import pdb; pdb.set_trace()
                    with open(pdbname+"_chain"+key+"_"+str(dummy)+"_"+key2+".pdb",'w') as f:
                        for line in pdb_by_chain[key]:
                            f.write(line)
                        for line in hetatoms[key2]:
                            f.write(line)

        # remove all the old files
        os.remove(pdbname+".pdb")


    def get_pdbs_from_PDB(self):
        '''
        downloads the pdbs from the PDB based on SMILES on ligands
        @return:
        '''

        url = 'http://www.rcsb.org/pdb/rest/search'

        queryText = """
<?xml version="1.0" encoding="UTF-8"?>

<orgPdbCompositeQuery version="1.0">

<queryRefinement>
<queryRefinementLevel>0</queryRefinementLevel> 

<orgPdbQuery>

<queryType>org.pdb.query.simple.ResolutionQuery</queryType>

<description>ResolutionQuery: refine.ls_d_res_high.comparator=between refine.ls_d_res_high.min=0.0 refine.ls_d_res_high.max="""+str(self.maxlength)+""" </description>

<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>

<refine.ls_d_res_high.min>0.0</refine.ls_d_res_high.min>

<refine.ls_d_res_high.max>"""+str(self.maxlength)+"""</refine.ls_d_res_high.max>

</orgPdbQuery> 

</queryRefinement>
  
<queryRefinement>                                                                                               
 
<queryRefinementLevel>1</queryRefinementLevel>  
<conjunctionType>and</conjunctionType> 
<orgPdbQuery>
<queryType>org.pdb.query.simple.ChemSmilesQuery</queryType>

<description>Chemical structure (SMILES/SMARTS): Structure is  and Smiles is """+self.smiles+""" and Search type is Similar and Similarity is 0.7 and Polymeric type is Any</description>

<structure></structure>

<smiles>"""+self.smiles+"""</smiles>

<searchType>Substructure</searchType>

<similarity>"""+str(self.similarity)+"""</similarity>

 <polymericType>Any</polymericType>

</orgPdbQuery>
</queryRefinement>

<queryRefinement>                                                                                                 
<queryRefinementLevel>2</queryRefinementLevel>                                                                    
<conjunctionType>and</conjunctionType>         

<orgPdbQuery>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
    <containsProtein>Y</containsProtein>
    <containsDna>N</containsDna>
    <containsRna>N</containsRna>
    <containsHybrid>N</containsHybrid>
  </orgPdbQuery>
</queryRefinement>

<orgPdbQuery>

<queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>

<description>Sequence Length Search : Min Sequence Length=4 Max Sequence Length="""+str(self.maxlength)+"""</description>

<v_sequence.chainLength.min>4</v_sequence.chainLength.min>

<v_sequence.chainLength.max>"""+str(self.maxlength)+"""</v_sequence.chainLength.max>

</orgPdbQuery>

</orgPdbCompositeQuery>

"""

        print "querying PDB...\n"
        req = urllib2.Request(url, data=queryText)

        f = urllib2.urlopen(req)

        # Contains pdb ids
        result = f.read().rstrip().split('\n')
        if result:

            print "Found number of PDB entries:", result.count('\n')

        else:

            print "Failed to retrieve results"


        for i in result:
            tmp = self.convert_to_pdb_gz(i)
            print "Downloading PDB structure '%s'..." % i
            self.get_gz_pdbfile(i)

    def convert_pdb_smi(self):
        # ~greisen/ExternalProgram/openbabel_selfcompiled/bin/babel -ipdb hcy.pdb -osmi query.smi
        ###import pdb; pdb.set_trace()
        exe = self.babelbin+self.babel+" -i"+self.format+" "+self.fragment_file+" -osmi query.smi"
        print exe

        subprocess.Popen(exe,shell=True).wait()

    def set_smiles(self):
        with open("query.smi") as f:
            for line in f:
                self.smiles = line.split()[0]
                self.smiles.strip()

        print "The PDB will be searched with this chemical fragment: "+self.smiles+"\n"

    def align_pdb_obfit(self):
        init_pdb = ""
        align_to_this = 0
        pdbs = os.listdir("./")
        for pdb in pdbs:
            if(pdb.endswith(".pdb") and align_to_this == 0):
                init_pdb = pdb
                align_to_this = 1
            elif(pdb.endswith(".pdb")):
                exe = self.obfit+" \""+self.smiles+"\" "+self.fragment_file+" "+pdb+"> "+pdb.split()[0]+"_aligned.pdb"
                p = subprocess.Popen(exe,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                out,output = p.communicate()
                if( len(output) == 0 and pdb != self.fragment_file):
                    exe_rm_1 = "rm -f "+pdb
                    exe_rm_2 = "rm -f "+pdb.split()[0]+"_aligned.pdb"

                    subprocess.Popen(exe_rm_1,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    subprocess.Popen(exe_rm_2,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

                '''
                # test 28-05-2015
                # exe = self.obfit+" \""+self.smiles+"\" "+init_pdb+" "+pdb+"> "+pdb.split()[0]+"_aligned.pdb"
                exe = self.obfit+" \""+self.smiles+"\" "+self.fragment_file+" "+pdb+"> "+pdb.split()[0]+"_aligned.pdb"
                ## print exe
                subprocess.Popen(exe,shell=True).wait()
                '''


            else:
                continue




    def main(self):

        parser = argparse.ArgumentParser(description="Takes a pdb with a chemical fragment and seraches the PDB for the occurance of this fragment.")
        # get the initial rosetta design as input
        parser.add_argument("-s", dest="fragment_file", help="This file contains the coordinates of the chemical fragment" )
        parser.add_argument("--maxlength", dest="maxlength", help="The max length of the protein (Default=200 aa)", default="200", type=str )
        parser.add_argument("--resolution", dest="resolution", help="Resolution of crystal structure (Default=2.0)", type=str, default="2.0" )
        parser.add_argument("--similarity", dest="similarity", help="The chemical similarity between the fragment searched (Default=1.0)", type=str )
        #parser.add_argument("--babel", dest="babel", help="The path to the executable for openbabel (Default - the dig system in the Bakerlab" )
        parser.add_argument("--format", dest="format", help="Format to convert from (Default=pdb)",default="pdb" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # converts the query pdb to smile format
        self.convert_pdb_smi()

        self.set_smiles()

        self.get_pdbs_from_PDB()

        self.align_pdb_obfit()


if __name__ == "__main__":
    run = SearchPDBwSmiles()
    run.main()
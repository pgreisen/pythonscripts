import urllib2
import urllib 
import gzip
import os
'''
How is this executed


'''


#smiles = '''C1CC[C@@H](C)[C@H](CCC)[C@@H]1C'''
# MPA 
##smiles = '''C1(=O)OCc2c(C)c(OC)cc(O)c12'''
##smiles = '''O(C)c1c(C)ccc(O)c1'''
# methoxy
#smiles='COC1=CC=CC=C1'
# HCY 
smiles='O[C@@H](C)C(=O)CO'

# Requires ent.gz file
# Returns pdb file and remove ent.gz file

# The chemical similarity between query and 
# target
similarity="0.99"

# maximum length of the protein
maxlength ="200"

# maximum resolution of the protein structure
maxresolution = "2.0"

def writes_to_pdb_format(filename,pdbname):
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
def convert_to_pdb_gz(pdbid):
    pdbid = pdbid.lower()
    filename = "pdb"+pdbid+".ent.gz"  

    return filename



# Requires Protein Data Bank url
# Writes pdb file
def get_gz_pdbfile(pdbname):
    gz_name = convert_to_pdb_gz(pdbname)

    url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'+gz_name[4:6]+'/'+gz_name

    urllib.urlretrieve(url, gz_name) 
    
    writes_to_pdb_format(gz_name,pdbname)



url = 'http://www.rcsb.org/pdb/rest/search'

queryText = """

<?xml version="1.0" encoding="UTF-8"?>


<orgPdbCompositeQuery version="1.0">  


<queryRefinement>                                                                                                
<queryRefinementLevel>0</queryRefinementLevel> 

<orgPdbQuery>

<queryType>org.pdb.query.simple.ResolutionQuery</queryType>

<description>ResolutionQuery: refine.ls_d_res_high.comparator=between refine.ls_d_res_high.min=0.0 refine.ls_d_res_high.max="""+str(maxlength)+""" </description>

<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>

<refine.ls_d_res_high.min>0.0</refine.ls_d_res_high.min>

<refine.ls_d_res_high.max>"""+str(maxlength)+"""</refine.ls_d_res_high.max>

</orgPdbQuery> 

</queryRefinement>
  
<queryRefinement>                                                                                               
 
<queryRefinementLevel>1</queryRefinementLevel>  
<conjunctionType>and</conjunctionType> 
<orgPdbQuery>
<queryType>org.pdb.query.simple.ChemSmilesQuery</queryType>

<description>Chemical structure (SMILES/SMARTS): Structure is  and Smiles is """+smiles+""" and Search type is Similar and Similarity is 0.7 and Polymeric type is Any</description>

<structure></structure>

<smiles>"""+smiles+"""</smiles>

<searchType>Substructure</searchType>

<similarity>"""+str(similarity)+"""</similarity>

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

<description>Sequence Length Search : Min Sequence Length=4 Max Sequence Length="""+str(maxlength)+"""</description>

<v_sequence.chainLength.min>4</v_sequence.chainLength.min>

<v_sequence.chainLength.max>"""+str(maxlength)+"""</v_sequence.chainLength.max>

</orgPdbQuery>

</orgPdbCompositeQuery>

"""


# print "query:\n", queryText

print "querying PDB...\n"

req = urllib2.Request(url, data=queryText)

f = urllib2.urlopen(req)

# Contains pdb ids
result = f.read().rstrip().split('\n')


if result:

    print "Found number of PDB entries:", result.count('\n')

else:

    print "Failed to retrieve results"


# correct url = "ftp://ftp.wwpdb.org"


# adding
# '/pub/pdb/data/structures/

# filename = "4HLW"
#filename = "pdb4hlw.ent.gz"
#pdb_dir = "divided" if not obsolete else "obsolete" 
#print pdb_dir

for i in result:
    tmp = convert_to_pdb_gz(i)
    print tmp



    print "Downloading PDB structure '%s'..." % i
    get_gz_pdbfile(i)
    #    get_gz_pdbfile(url,filename)


#urllib.urlretrieve(url, filename) 

#writes_to_pdb_format(filename)


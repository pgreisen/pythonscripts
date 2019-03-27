# Requires Protein Data Bank url
# Writes pdb file
import sys
pdbid=sys.argv[1]
from urllib.request import urlopen
import urllib.request
with urllib.request.urlopen("https://files.rcsb.org/download/"+pdbid.upper()+".pdb") as response:
   html = response.read()

a=html.decode('ascii')
b=a.split('\n')

with open(pdbid+".pdb",'w') as f:
   for i in b:
      f.write(i+"\n")
   

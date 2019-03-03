# Requires Protein Data Bank url
# Writes pdb file
org="EXPRESSION_SYSTEM:"
from urllib.request import urlopen
import urllib.request
with urllib.request.urlopen("https://files.rcsb.org/header/1QLU.pdb") as response:
   html = response.read()

#line = html.decode("utf-8")
#print(line)

a=html.decode('ascii')#.strip().split(',')
#print(a,type(a))
b=a.split('\n')
#b=list(a)
for i in b:
   print(i)
   if(org in i):
      print(i.split(':')[1])
#encoding="utf-8"
#output=html.decode()
#output=str(html, encoding)
#with open(html) as f:
#for ln in line:
#   print(ln)

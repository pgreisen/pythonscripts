import urllib.request
code = "P81277"
data = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + code + ".txt").read()
org="SQ"
a=data.decode('ascii')#.strip().split(',')
#print(a,type(a))
b=a.split('\n')
fastaseq=""
seq=False
organism=""
for i in b:
   if("OS" in i):
      organism=i.replace("OS   ","")
   if(seq):
      fastaseq += i.strip()
   if(org in i):
      seq=True
fastaseq = fastaseq.replace("/","")
fastaseq = fastaseq.replace(" ","")
print(code,organism)
print(fastaseq)



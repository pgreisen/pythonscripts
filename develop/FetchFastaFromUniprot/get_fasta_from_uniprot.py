import urllib.request

def get_fasta_from_uniprot(uniprot_id):
    data = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + uniprot_id + ".txt").read()
    org = "SQ"
    prt_name="DE   RecName: Full="
    a = data.decode('ascii')
    b = a.split('\n')
    fastaseq = ""
    seq = False
    organism = ""
    for i in b:
        if ("OS" in i):
            organism = i.replace("OS   ", "")
        if (seq):
            fastaseq += i.strip()
        if (org in i):
            seq = True
        if(prt_name in i):
            fastafunction = i.replace(prt_name,"")
    fastaseq = fastaseq.replace("/", "")
    fastaseq = fastaseq.replace(" ", "")
    fastaid=code+"_"+organism
    return fastaid, fastaseq, fastafunction

code = "P81277"
fastaid, fastaseq, fastafunction = get_fasta_from_uniprot(code)
print(fastaid)
print(fastaseq)
print(fastafunction)

import sys

def get_seq(fastafile):
    seq = ""
    with open(fastafile,'r') as f:
        for line in f:
            if(line[0] != ">" ):
                seq += line.strip()
    return seq
                
fastafile = sys.argv[1]
aa_list = ['A','C','H','D','E','K','S','R','W','N','Q','V','I','L','F','M','T','P','G','Y']
header="total 19\n"
seq = get_seq( fastafile )

for pos in range(len(seq) ):
    name = seq[pos]+"_"+str(pos)+".mutfile"
    wt_aa = seq[pos]
    with open( name,'w' ) as f:
        f.write(header)
        for i in aa_list:
            if( i != seq[pos] ):
                f.write("1\n"+seq[pos] +"  "+str(pos+1)+"  "+i+"\n")

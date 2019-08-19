# get fasta file into dictionary
FASTAFILE = "sciworm.fasta"
seqs = {}
with open(FASTAFILE, 'r') as f:
    for i in f:
        if(i[0] == ">"):
            key = i.strip()[1:]
            seqs[key] = ""
        else:
            seqs[key] += i.strip()
df = pd.DataFrame.from_dict(seqs,orient='index').reset_index().rename(columns={'index' : 'Name', 0 : 'sequence'})

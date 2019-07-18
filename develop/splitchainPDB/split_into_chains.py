import sys,os


chain_breaker = ["TER"]
chains = ['A','B','C','D','E','F']

def set_chains(pdbfile,dictionary):
    key = ""
    pdbname = pdbfile.split(".")[0]
    

    with open(pdbfile) as f:
        for line in f:
            if(line[0:4] == "ATOM" or line[0:4] == "HETA" ):
                chainid_ = line[21:22].strip()
                key = pdbname+"_chain_"+chainid_
                if( key not in dictionary.keys() ):
                    dictionary[key] = []
                dictionary[key].append(line)

sequences = {}
pdbfile = sys.argv[1]
set_chains( pdbfile, sequences )

# print sequences
for key in sequences:
    with open(key+".pdb", 'w') as f:
        for line in sequences[key]:
            f.write( line )

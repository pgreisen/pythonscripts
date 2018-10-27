import sys,os


chain_breaker = ["TER"]

def set_chains(pdbfile,dictionary):

    with open(pdbfile) as f:
        key = ""

        dummy = 0
        for line in f:
            if (line[0:3] in chain_breaker or dummy == 0):
                key = pdbfile[0:4].replace('.','') + "_" + str(dummy) + "_chain"
                dummy += 1


            if(line[0:4] == "ATOM" or line[0:4] == "HETA" ):
                if(key not in dictionary.keys() ):
                    dictionary[key] = []
                dictionary[key].append(line)



sequences = {}
pdbfile = sys.argv[1]
set_chains( pdbfile, sequences )


pdbname = pdbfile.split(".")[0]
# print sequences
for key in sequences:
    with open(key+".pdb", 'w') as f:
        for line in sequences[key]:
            f.write( line )

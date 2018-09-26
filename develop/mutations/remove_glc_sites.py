import sys,re

datafile = sys.argv[1]
sequences = {}

with open(datafile, 'r') as f:
    for line in f:

        if (line[0] == ">"):
            tmp = line.split()
            ids = tmp[0]
            sequences[ids] = ""
        else:
            tmp = line.strip()
            sequences[ids] = sequences[ids] + tmp


with open("glc_removed_"+datafile, 'w') as f:
    remove_ = 0

    for key in sequences:
        # test if sequence contain glc site
        result = re.search(r'N[A-Z][ST]', sequences[key] )
        chemical_stability = re.search(r'[ND][GAS]', sequences[key] )
        if(result is not None):
            remove_ += 1
            continue
        #elif(chemical_stability is not None):
        #    remove_ += 1
        else:
            f.write(key + '\n')
            f.write(sequences[key] + '\n')

print "Number of sequences removed: ",remove_





import numpy as np

def get_min_max_gene_length(fastafile):
    var_min_max = {}
    with open(fastafile, 'r') as f:
        for line in f:
            if(line[0] == '>'):
                tmpline = line[1:].split("_")
                for i in tmpline[:-1]:
                    min_max_.append(int(i)) 
                var_min_max[line[1:] = (min(min_max_),max(min_max_))
    return var_min_max


# the length of the synthesized gene
gene_construct_length = 100

gene_length = 1500

# hash of designs
# collect beginning and end of each design
'''
a = {
var1 : (min,max)
}


'''
# Construct clusters

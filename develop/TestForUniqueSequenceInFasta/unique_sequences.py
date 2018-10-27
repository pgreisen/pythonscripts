import sys

'''
Test if all fasta sequences in a file are unique and print the header for 
those which are not unique.

Script is executed as 

python unique_sequences.py FASTAFILE




'''

def get_fastafiles(datafile):
    sequences = {}
    dummy = 0
    with open(datafile, 'r') as f:
        for line in f:
            if (line[0] == ">"):
                tmp = line.split()
                ids = tmp[0]+"#"+str(dummy)
                sequences[ids] = ""
                dummy += 1
            else:
                tmp = line.strip()
                sequences[ids] = sequences[ids] + tmp
    print "Total number of fasta sequences: ", len(sequences.keys())
    return sequences


def get_unique(sequences):
    uniq_seq = []
    for key in sequences:
        uniq_seq.append(sequences[key])
    print "Total number of fasta sequences: ", len(uniq_seq)
    print "Unique number of fasta sequences: ", len(set(uniq_seq))
    return uniq_seq

def get_doublicates(sequences):
    reverse_dic = {}
    for k, v in sequences.iteritems():
        reverse_dic[v] = reverse_dic.get(v, [])
        reverse_dic[v].append(k)

    for k, v in reverse_dic.iteritems():
        if (len(reverse_dic[k]) > 1):
            print v
            for i in v[1:]:
                sequences.pop(i, None)
    with open("unique.fasta",'w') as f:
        for k,v in sequences.iteritems():
            import re
            k =  re.sub(r'\#.*$', '', k)
            f.write(k+'\n')
            f.write(v+'\n')
    return sequences


datafile = sys.argv[1]

fasta_ = get_fastafiles(datafile)

uniq_seq = get_unique(fasta_)

if(len(uniq_seq) != len(set(uniq_seq))):
    print "Removing doublicates from file"
    reverse_dic = get_doublicates(fasta_)
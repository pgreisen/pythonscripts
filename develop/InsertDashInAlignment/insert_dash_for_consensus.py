import sys

def write_to_file(filename, alignment_with_x):
    # number of sequences in alignments
    alignment_size = len(alignment_with_x[0])
    # write to file
    with open(filename, 'w') as f:
        for alignment_nr in range(alignment_size):
            f.write(">" + str(alignment_nr) + "\n")
            for key in alignment_with_x.keys():
                f.write(alignment_with_x[key][alignment_nr])
            f.write("\n")


datafile = sys.argv[1]
sequences = {}
pos_sequences = {}


with open(datafile,'r') as f:
    for line in f:
        if(line[0] == ">"):
            tmp = line.split()
            ids = tmp[0]
            sequences[ids] = ""
        else:
            tmp = line.strip()
            sequences[ids] = sequences[ids] + tmp


alignment_with_x_pos = {}

prt_length = 0

for i in sequences.keys():
    prt = list(sequences[i])
    prt_length = len(sequences[i])
    for aa in range(prt_length):
        if (aa not in alignment_with_x_pos.keys()):
            alignment_with_x_pos[aa] = []
        alignment_with_x_pos[aa].append(prt[aa])


for tmp in alignment_with_x_pos.keys():
    if (len(set(alignment_with_x_pos[tmp])) == 1):
            alignment_with_x_pos[tmp] = [aa_sub.replace(aa_sub, "-") for aa_sub in alignment_with_x_pos[tmp]]

write_to_file("LC_positive.fasta", alignment_with_x_pos)








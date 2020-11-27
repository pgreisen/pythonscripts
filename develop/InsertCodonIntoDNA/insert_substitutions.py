import sys

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


aa_codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}

insert_mutations ={
    '130' : 'G',
    '131' : 'G',
    '132' : 'G',
    '133' : 'G',
    '134' : 'G',
    '135' : 'G',
    '136' : 'G',
    '137' : 'G',
    '138' : 'G',
    '139' : 'G',
}


# a function to translate a single codon
def translate_codon(codon):
    return gencode.get(codon.upper(), 'x')
 
# a function to split a sequence into codons
def split_into_codons(dna, frame):
    codons = []
    for i in range(frame - 1, len(dna)-2, 3):
        codon = dna[i:i+3]
        codons.append(codon)
    return codons
 
# a function to translate a dna sequence in a single frame
def translate_dna_single(dna, frame=1):
    codons = split_into_codons(dna, frame)
    amino_acids = ''
    for codon in codons:
        amino_acids = amino_acids + translate_codon(codon)
    return amino_acids
 
# a function to translate a dna sequence in 3 forward frames
def translate_dna(dna):
    all_translations = []
    for frame in range(1,4):
        all_translations.append(translate_dna_single(dna, frame))
    return all_translations

def get_fasta_design_seq(fastafile):
    ff = {}
    with open(fastafile,'r') as f:
        for line in f:
            if(line[0] == '>'):
                key = line.split()[0]
                ff[key] = ""
            else:
                ff[key]= ff[key] + line.strip()
    return ff

def insert_new_codon(hash_of_codons,chains,chainid=0):
    '''
    hash_of_codons : dictonary with position as key and codon as value
    chains : dictionary with fasta files
    key : sequence to work on where key is identifier in hash

    '''
    newchain = chains[chainid]
    for key in hash_of_codons:
        # -1 due to newchain being a list
        newchain[int(key)-1] = hash_of_codons[key]
    return newchain

def get_codons_and_positions(insertions):
    aa_codons = {}
    with open(insertions) as f:
        for line in f:
            pos_,codon_ = line.strip().split(',')
            aa_codons[pos_] = codon_
    return aa_codons


def write2file(newchain,name):
    with open(name+".fasta",'w') as f:
        f.write(">"+name+"\n")
        for i in newchain:
            f.write(i)

fastafile = sys.argv[1]

# A text file with position,codon
insertions = sys.argv[2]

new_codons_to_insert = get_codons_and_positions(insertions)
print(new_codons_to_insert)


# return an dictionary with key and sequence
dna_seq = get_fasta_design_seq(fastafile)

chains = {}
for key in dna_seq:
    # return a list with codons from fasta file
    chains[key] = split_into_codons(dna_seq[key], 1)

new_chains = {}

for key in chains.keys():
    new_chains[key] = insert_new_codon(new_codons_to_insert,chains,key)

for i in new_chains.keys():
    write2file(new_chains[i],i[1:]+"_library")

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


'''


'''

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


pte = "ggaaggatttcagaattcatcaccaacagcggcgaccgtatcaacaccgtccgtggtccgatcaccatctctgaggcgggcttcaccctgactcacgaacacatttgcggttctagcgcaggttttctgcgcgcttggccggagtttttcggttctcgtgctgctctggtggaaaaagcggttcgtggcctgcgtcgtgcgcgtgcggctggtgtgcgtaccatcgttgacgtttctaccttcgacattggtcgtgatgtttctctgctggccgaggtttctgaagcggccgatgttcacattgttgcagcgactggtctgtgggaagatccgccgctgtctatgcgtctgcgctctgttgaagaactcacccagttctttctccgtgaaatccagtacggcatcgaggacacgggtatccgtgccggtatcattaaagttgccaccaacggtaaagcgaccccgtttcaggaactggttctgcgtgcagcagctcgtgcctccctcgccaccggcgttccggtcaccacccacaccgacgcttctcagcgtgacggtgaacagcaggcggcgatcttcgaaagcgaaggtctggacccgtctcgtgtttgtatcggtcactctgacgacaccgatgatctggactacctgaccgcgctcgcggctcgtggttacctgattggcctggatggtattccgcactctgcgatcggcctcgaagacaacgcatctgcgtccgctctgctcggtaatcgctcttggcagacccgtgcgctgctgatcaaagcgctgatcgaccagggctacgttaaacagatcctggtttctaacgattggctgttcggtttttcttcttgggttaccaacatcatggacgttatggactctgttaacccagacggtatggcgttcatcccgctgcgtgttatcccgttcctgcgcgagaaaggtgttccacaagagacgctggcgaccatcaccgttgaaaaccctgctcgtttcctgtctccaaccctccgtgcttcttgataactgcaggcaagcttggc"


# GRISEFITNSG
# The protein is kept between these sequencs
# RAS

test = translate_dna(pte)
#print test
print test[0][11:338]

##print pte[33:1014]
pte_dna = pte[33:1014]
test2 = translate_dna(pte[33:1014])
##print test2[0], len(test2[0])
##print "Resi 132 is ",test2[0][96:98]
##print "AA: ",aa_codon_table['W'][0]
# insert the mutations into the string
new_chain = pte_dna
for key in insert_mutations:
    # key is residue number
    # multiple by 3 to get sequence
    # print key,aa_codon_table[ insert_mutations[key]  ][0]

    insert_here = 3*(int( key) - 35)
    print insert_here
    new_chain = new_chain[0:insert_here]+aa_codon_table[ insert_mutations[key]  ][0]+new_chain[insert_here+3:]

    #print insert_mutations[key], key, aa_codon_table[ insert_mutations[key]  ][0]
    #print test2[0][0:int(key)]
    #print aa_codon_table[ insert_mutations[key]  ][0]
    # print translate_dna( aa_codon_table[ insert_mutations[key]  ][0] )[0]

print translate_dna( new_chain )[0]

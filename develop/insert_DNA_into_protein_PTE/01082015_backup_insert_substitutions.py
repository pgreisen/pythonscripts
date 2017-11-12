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
    ff = ""
    with open(fastafile,'r') as f:
        for line in f:
            if(line[0] != '>'):
                print line
                ff += line.strip()
    return ff



fastafile = sys.argv[1]
ff = get_fasta_design_seq(fastafile)

pte = "ggaaggatttcagaattcatcaccaacagcggcgaccgtatcaacaccgtccgtggtccgatcaccatctctgaggcgggcttcaccctgactcacgaacacatttgcggttctagcgcaggttttctgcgcgcttggccggagtttttcggttctcgtgctgctctggtggaaaaagcggttcgtggcctgcgtcgtgcgcgtgcggctggtgtgcgtaccatcgttgacgtttctaccttcgacattggtcgtgatgtttctctgctggccgaggtttctgaagcggccgatgttcacattgttgcagcgactggtctgtgggaagatccgccgctgtctatgcgtctgcgctctgttgaagaactcacccagttctttctccgtgaaatccagtacggcatcgaggacacgggtatccgtgccggtatcattaaagttgccaccaacggtaaagcgaccccgtttcaggaactggttctgcgtgcagcagctcgtgcctccctcgccaccggcgttccggtcaccacccacaccgacgcttctcagcgtgacggtgaacagcaggcggcgatcttcgaaagcgaaggtctggacccgtctcgtgtttgtatcggtcactctgacgacaccgatgatctggactacctgaccgcgctcgcggctcgtggttacctgattggcctggatggtattccgcactctgcgatcggcctcgaagacaacgcatctgcgtccgctctgctcggtaatcgctcttggcagacccgtgcgctgctgatcaaagcgctgatcgaccagggctacgttaaacagatcctggtttctaacgattggctgttcggtttttcttcttgggttaccaacatcatggacgttatggactctgttaacccagacggtatggcgttcatcccgctgcgtgttatcccgttcctgcgcgagaaaggtgttccacaagagacgctggcgaccatcaccgttgaaaaccctgctcgtttcctgtctccaaccctccgtgcttcttgataactgcaggcaagcttggc"

##test = translate_dna(pte)
#print test
##print test[0][11:338]

##print pte[33:1014]
pte_dna = pte[33:1017]
test2 = translate_dna(pte[33:1017])

print pte_dna

# insert the mutations into the string
new_chain = pte_dna
for key in insert_mutations:
    # key is residue number
    # multiple by 3 to get sequence
    # print key,aa_codon_table[ insert_mutations[key]  ][0]
    insert_here = 3*(int( key) - 35)
    ###print "insertion is performed here:",insert_here, key
    new_chain = new_chain[0:insert_here]+aa_codon_table[ insert_mutations[key]  ][0]+new_chain[insert_here+3:]

moshe_gene = translate_dna( new_chain )[0]
#print "MOSHE: ",moshe_gene, len( moshe_gene )
#print "DESIGN: ",ff ,len( ff )
assert len( moshe_gene ) == len( ff )

design_pte_seq = pte_dna
designed_fasta_seq = ""
for aa in range ( len( moshe_gene ) ):
    if( moshe_gene[aa] != ff[aa] ):
        end = (aa + 1)*3
        begin = end - 3
        dna_sub = aa_codon_table[ff[aa] ][0].lower()
        design_pte_seq = design_pte_seq[0:begin]+dna_sub+design_pte_seq[end:]
    else:
        designed_fasta_seq += moshe_gene[aa]


##print translate_dna( design_pte_seq )[0]
final_seq = pte[0:33]+design_pte_seq+pte[1018:]
##print len( final_seq )
print final_seq

dna_file = fastafile.split('.')[0]

with open(dna_file+'.gb','w') as f:
    for line in final_seq:
        f.write(line)

import sys,argparse, random

class InsertSubstitutions:

    def __init__(self):
        self.DNAseq = ""
        self.DNAfile = ""
        self.codon_table = {}
        self.codon_table_file = ""
        self.mutations = ""
        self.mutational_table = {}

        self.aa = {
            'H': 'HIS',
            'E': 'GLU',
            'D': 'ASP',
            'C': 'CYS',
            'K': 'LYS',
            'T': 'THR',
            'S': 'SER',
            'R': 'ARG',
            'Q': 'GLN',
            'N': 'ASN',
            'Y': 'TYR',
            'M': 'MET',
            'P': 'PRO',
            'A': 'ALA',
            'V': 'VAL',
            'L': 'LEU',
            'I': 'ILE',
            'F': 'PHE',
            'W': 'TRP',
            'G': 'GLY'
        }


    def setup_codon_table(self):

        with open(self.codon_table_file,'r') as f:
            for line in f:
                tmpline = line.split()
                tmpaa = tmpline[0].upper()
                if( tmpaa not in self.codon_table):
                    self.codon_table[tmpaa] = {}
                    self.codon_table[tmpaa][tmpline[1]] = round(float(tmpline[2]),2)
                else:
                    self.codon_table[tmpaa][tmpline[1]] = round(float(tmpline[2]),2)


    def get_codon(self,aa):
        lngth = len(self.codon_table[self.aa[aa]])
        rnr = random.random()
        sum = 0
        for key, value in self.codon_table[self.aa[aa]].items():
            sum += value
            if(rnr < sum):
                return key
        return key


    # a function to translate a single codon
    def translate_codon(self, codon):
        return gencode.get(codon.upper(), 'x')

    # a function to split a sequence into codons
    def split_into_codons(self, dna, frame):
        codons = []
        for i in range(frame - 1, len(dna) - 2, 3):
            codon = dna[i:i + 3]
            codons.append(codon)
        return codons

    # a function to translate a dna sequence in a single frame
    def translate_dna_single(self, dna, frame=1):
        codons = split_into_codons(dna, frame)
        amino_acids = ''
        for codon in codons:
            amino_acids = amino_acids + translate_codon(codon)
        return amino_acids

    # a function to translate a dna sequence in 3 forward frames
    def translate_dna(self, dna):
        all_translations = []
        for frame in range(1, 4):
            all_translations.append(translate_dna_single(dna, frame))
        return all_translations

    def get_fasta_design_seq(self, fastafile):
        ff = ""
        with open(fastafile, 'r') as f:
            for line in f:
                if (line[0] != '>'):

                    ff += line.strip()
        return ff

    def insert_mutation(self, DNAseq, mutational_dictionary):
        new_chain = DNAseq


        for key in mutational_dictionary.keys():
            insert_here = 3 * (int(key) -1)
            new_chain = new_chain[0:insert_here] + self.get_codon(mutational_dictionary[key]).lower() + new_chain[insert_here + 3:]
        return new_chain


    def set_mutational_table(self):
        tmpstring = self.mutations.split('_')
        for i in tmpstring:
            # position
            key = i[1:-1]
            # amino acids
            mut = i[-1]
            self.mutational_table[key] = mut



    def get_mutational_table(self,mutations):
        tmpstring = mutations.split('_')
        mutational_table = {}
        for i in tmpstring:
            # position
            key = i[1:-1]
            # amino acids
            mut = i[-1]
            mutational_table[key] = mut
        return mutational_table


    def to_file(self,fastaheader,sequence):
        with open(fastaheader+".fasta",'w') as f:
            f.write(">"+fastaheader+"\n")
            f.write(sequence)

    def run(self, DNAfile, codonFILE, mutations):

        DNAseq = self.get_fasta_design_seq(DNAfile)
        self.codon_table_file = codonFILE
        self.mutations = mutations

        self.setup_codon_table()
        mutational_table = self.get_mutational_table(mutations)
        newDNA = self.insert_mutation(DNAseq, mutational_table)

        # Make sure that no deletions have been made
        assert len(newDNA) == len(DNAseq)

        return self.mutations, newDNA



    def main(self):
        parser = argparse.ArgumentParser(description="")
        # get the initial rosetta design as input
        parser.add_argument("--dna", dest="DNAfile", help="File with template DNA")
        parser.add_argument("--codon", dest="codon_table_file", help="Codon scoring table")
        parser.add_argument("--mutations", dest="mutations", help="Mutational string of the form A1X_B2Y etc")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.DNAseq = self.get_fasta_design_seq(self.DNAfile)
        self.setup_codon_table()

        self.set_mutational_table()
        newDNA = self.insert_mutation(self.mutational_table)

        # Make sure that no deletions have been made
        print(newDNA, DNAseq)
        assert len(newDNA) == len(self.DNAseq)

        self.to_file(self.mutations, newDNA)


if __name__ == "__main__":
   run = InsertSubstitutions()
   run.main()



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
        print(self.codon_table)

    def get_codon(self,aa):
        lngth = len(self.codon_table[aa])
        rnr = random.random()
        sum = 0
        for key, value in self.codon_table[aa].items():
            sum += value
            if(rnr < sum):
                return value
        return value


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
                    print(line)
                    ff += line.strip()
        return ff

    def insert_mutation(self, mutational_dictionary):
        new_chain = self.DNAseq
        for key in mutational_dictionary.keys():
            insert_here = 3 * (int(key) - 35)
            #new_chain = new_chain[0:insert_here] + aa_codon_table[insert_mutations[key]][0] + new_chain[insert_here + 3:]
            new_chain = new_chain[0:insert_here] + self.get_codon(mutational_dictionary[key]) + new_chain[insert_here + 3:]

        return new_chain


    def set_mutational_table(self):
        tmpstring = self.mutations.split('_')
        for i in tmpstring:
            # position
            key = i[1:-1]
            # amino acids
            mut = i[-1]
            self.mutational_table[key] = mut


    def main(self):
        parser = argparse.ArgumentParser(description="")
        # get the initial rosetta design as input
        parser.add_argument("--dna", dest="DNAfile", help="File with template DNA")
        parser.add_argument("--codon", dest="codon_table_file", help="Codon scoring table")
        parser.add_argument("--mutations", dest="mutations", help="Mutational string of the form A1X_B2Y etc")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.setup_codon_table()
        self.set_mutational_table()
        newDNA = self.insert_mutation(self.mutational_table)
        print newDNA

if __name__ == "__main__":
   run = InsertSubstitutions()
   run.main()



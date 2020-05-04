# Multiple site ddG
#
# Input pdb file 
# Output file 

import argparse

class GenerateMutationfileDdg:


    def __init__(self):
        self.pdbfile=""
        self.pdbfile_singleletter = ""
        self.modulus = 10
        self.pdb_ddg = {}
        self.aas = "ARNDCQEGHILKMFPSTWYV"

    def get_single_letter_AA(self,residue):
        aa = {
            "ALA": "A",
            "ILE": "I",
            "LEU": "L",
            "VAL": "V",
            "MET": "M",
            "PHE": "F",
            "TYR": "Y",
            "ARG": "R",
            "LYS": "K",
            "TRP": "W",
            "ASN": "N",
            "GLN": "Q",
            "ASP": "D",
            "GLU": "E",
            "SER": "S",
            "THR": "T",
            "GLY": "G",
            "PRO": "P",
            "CYS": "C",
            "HIS": "H"
        }
        if (isinstance(residue, str)):
            return aa[residue]
        return "NaN"

    def set_single_letter_pdb(self):
        with open(self.pdbfile,'r') as f:
            for line in f:
                if(line[0:4] == "ATOM"):
                    if( line[13:15] == "CA"):
                        self.pdbfile_singleletter += self.get_single_letter_AA(line[17:20].strip())
        print(self.pdbfile_singleletter)


    def generate_mutfile(self):
        for i in range(len(self.pdbfile_singleletter)):
            key=self.pdbfile_singleletter[i]+"_"+str(i+1)
            self.pdb_ddg[key] = i+1


    def write_mutation_file(self):
        # header = "total 1\n1\n"
        # 'D_178': 178
        total_count=0

        for key in self.pdb_ddg.keys():
            totsub = 0
            dummy = 0
            wt,pos = key.split("_")
            aas_ = self.aas.replace(wt,"")
            for aa_ in aas_:
                if(dummy % self.modulus == 0):
                    totsub += self.modulus
                    if(totsub > len(aas_)):
                        f.close()
                        tot_substitutions_per_file= len(aas_) - (totsub - self.modulus)
                    else:
                        tot_substitutions_per_file = self.modulus
                    f= open(wt + str(totsub) + "X_"+str(pos)+".mutfile", 'w')
                    header="total "+str(tot_substitutions_per_file)+"\n"
                    f.write(header)
                f.write("1\n")
                f.write(wt + " " + str(pos) + " " + aa_+"\n")
                dummy += 1


    def main(self):
        parser = argparse.ArgumentParser(description="Generate multiple files as input for ddG calculations in Rosetta")
        parser.add_argument('--pdbfile','-p', dest="pdbfile", help="pdbfile")
        parser.add_argument('--modulus','-m', dest="modulus", help="split each position into multiples")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
          setattr(self, item, args_dict[item])
        self.set_single_letter_pdb( )
        self.generate_mutfile()
        self.write_mutation_file()



if __name__ == "__main__":
    run = GenerateMutationfileDdg()
    run.main()






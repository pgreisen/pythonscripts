# Multiple site ddG
#
# Input pdb file 
# Output file 

import argparse

class GenerateMutationfileDdg:


    def __init__(self):
        self.pdbfile=""
        self.pdbfile_singleletter = ""
        self.modulus = 1
        self.pdb_ddg = {}
        self.aas = "ARNDCQEGHILKMFPSTWYV"
        self.mutations = ""
        self.muts = {}
        self.number = 1
        self.ssm = 0
        
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


    def generate_mutfile(self):
        # self.muts : list with mutations
        from itertools import combinations
        pos_ = self.muts.keys()
        combs = {}
        for i in range(1,self.number):
            key="N"+str(i)
            combs[key] = combinations(pos_, i)
        return combs
        


    def write_mutation_file(self, list_of_variants):
        # header = "total 1\n1\n"
        # 'D_178': 178
        total_count=0
        tmpsum = 0 
        sub = ""
        dummy = 0
        for tot_var in list_of_variants:
            totsub = 0
            dummy += 1
            nr = len(tot_var)
            tmpsum += nr
            sub += str(nr)+"\n"
            for var in tot_var:
                wt=var[0]
                pos=var[1:-1]
                aa=var[-1]
                sub += wt+" "+pos+" "+aa+"\n"
            #print(sub, tmpsum, self.modulus)
            if( tmpsum >= self.modulus  or dummy == len(list_of_variants) ):                
                f = open( str(total_count)+".mutfile", 'w')
                header="total "+str(tmpsum)+"\n"
                f.write(header)
                f.write(sub)
                f.close()
                total_count += 1
                tmpsum = 0
                sub = ""

##
    def write_mutation_file_ssm(self, list_of_variants):
        # header = "total 1\n1\n"
        # 'D_178': 178
        total_count=0
        tmpsum = 0 
        sub = ""
        dummy = 0
        nr = 1
        for tot_var in list_of_variants:
            totsub = 0
            dummy += 1
            tmpsum += nr
            sub += str(nr)+"\n"
            wt=tot_var[0]
            pos=tot_var[1:-1]
            aa=tot_var[-1]
            
            sub += wt+" "+pos+" "+aa+"\n"
            if( tmpsum >= self.modulus or dummy == len(list_of_variants) ):                
                f = open( str(total_count)+".mutfile", 'w')
                header="total "+str(tmpsum)+"\n"
                f.write(header)
                f.write(sub)
                f.close()
                total_count += 1
                tmpsum = 0
                sub = ""


                
    def generate_ssm(self):
        print(len(self.pdbfile_singleletter))
        for i in range(len(self.pdbfile_singleletter)):
            key=i+1
            if(key not in self.muts.keys()):
                self.muts[key] = []
            wt = self.pdbfile_singleletter[i]

            for aa in self.aas:
                if(aa == wt):
                    continue
                tmpline = wt+str(key)+aa
                self.muts[key].append(tmpline)

    def set_mutations(self):
        with open(self.mutations,'r') as f:
            for line in f:
                tmpline = line.strip()
                key = tmpline[1:-1]
                if(key not in self.muts.keys()):
                    self.muts[key] = []
                self.muts[key].append(tmpline)

    def main(self):
        parser = argparse.ArgumentParser(description="Generate multiple files as input for ddG calculations in Rosetta")
        parser.add_argument('--pdbfile','-p', dest="pdbfile", help="pdbfile")
        parser.add_argument('--modulus','-m', dest="modulus", help="split each position into multiples", default=1,type=int)
        parser.add_argument('--file','-f', dest="mutations", help="File with mutations e.g. A203E\nD325F etc")
        parser.add_argument('--number','-n', dest="number", help="Number of mutations e.g. single and double mutations with n = 2",default=1,type=int)
        parser.add_argument('--ssm','-x', dest="ssm", help="Perform full in silico saturation - default is false (0)",default=0,type=int)
        
        from itertools import combinations,product
        args_dict = vars(parser.parse_args())
        for item in args_dict:
          setattr(self, item, args_dict[item])

        list_of_variants = []
        if( self.ssm != 0  ):
            self.set_single_letter_pdb()
            self.generate_ssm()

            for i in self.muts.keys():
                print(self.muts[i])
                for j in self.muts[i]:
                    list_of_variants.append( j )
            self.write_mutation_file_ssm( list_of_variants)

        else:
            self.set_mutations()
            list = self.generate_mutfile()


            for i in list:
                for j in list[i]:

                    '''
                    list 1 N2 ('325', '23')
                    325 ['V325G']
                    23 ['H23W']
                    '''
                    b = []
                    for k in j:
                        b.append( self.muts[k] )
                    c = product(*b)
                    for l in c:
                        list_of_variants.append(l)
            self.write_mutation_file( list_of_variants)

if __name__ == "__main__":
    run = GenerateMutationfileDdg()
    run.main()

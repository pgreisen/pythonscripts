import sys, shutil, os, subprocess, argparse

'''

Generate input files for FoldX by splitting up the protein chain into multiple inputfiles

@require pdbfile

'''


class InputFileFoldX:

    def __init__(self):
        self.pdbfile = None
        self.protein_index = [] #set()
        self.protein_length = 0
        self.min = 0
        self.max = 0
        self.chain = 'A'
        self.fixed_residues = [""]
        # self.fixed_residues = ["55","57","169","201","230","301"]



    def get_string_of_mutations(self):
        print ""


    def get_single_residue(self,residuename):
        res_map = {
            "ALA" : 'A',
            "CYS" : 'C',
            "HIS" : 'H',
            "ASP" : 'D',
            "GLU" : 'E',
            "LYS" : 'K',
            "KCX" : 'K',
            "SER" : 'S',
            "ARG" : 'R',
            "TRP" : 'W',
            "ASN" : 'N',
            "GLN" : 'Q',
            "VAL" : 'V',
            "ILE" : 'I',
            "LEU" : 'L',
            "PHE" : 'F',
            "MET" : 'M',
            "THR" : 'T',
            "PRO" : 'P',
            "GLY" : 'G',
            "TYR" : 'Y'
            }
        return res_map[residuename]


#DA35a,RA36a,IA37a;

    def get_pdb_start_end(self):

        with open(self.pdbfile,'r') as f:
            for line in f:
                if ( line[0:4] == "ATOM"):
                    if ( line[13:15] == "CA"):
                        #chain = line[21:22]
                        #resname = self.get_single_residue( line[17:20])
                        residue_number = str(line[22:26]).strip()
                        if( residue_number not in self.fixed_residues ):
                            pdbstring = str(line[22:26]).strip()
                        else:
                            continue

                        self.protein_index.append( pdbstring )

        self.protein_length = len(self.protein_index)

        #self.min = min(self.protein_index)
        #self.max = max(self.protein_index)

    def get_template(self,mutations_string):

        template = '''NATRO
start
        '''+mutations_string
        return template



    def write_to_file( self, template, dummy ):
        with open(str(dummy)+"_resfile",'w') as f:
            for line in template:
                f.write( line)



    def main(self):

        parser = argparse.ArgumentParser(description=" Generates input file for ddG rosetta to compute mutational scan ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="pdbfile", help="Initial starting pdb")

        parser.add_argument("-p", "--processors", dest="processors", help="The number of processors to split the jobs into")

        parser.add_argument("-c", "--chain", dest="chain", help="Which chain to use", default='A')

        input_variables = parser.parse_args()

        self.pdbfile = input_variables.pdbfile

        self.chain = input_variables.chain

        self.get_pdb_start_end()

        dummy = 0

        split_into = self.protein_length / int(input_variables.processors)

        print "Length of protein is: ",self.protein_length
        last = False
        last_mutations_string = ""
        for i in range(1, self.protein_length+5, split_into):

            if( (len( self.protein_index )/split_into) > 0.0 and len(self.protein_index)> 4.9 ):
                mutations_string = self.protein_index.pop(0)+" "+self.chain+" ALLAA\n\t"+self.protein_index.pop(0)+" "+self.chain+" ALLAA\n\t"+self.protein_index.pop(0)+" "+self.chain+" ALLAA\n\t"+self.protein_index.pop(0)+" "+self.chain+" ALLAA\n\t"+self.protein_index.pop(0)+" "+self.chain+" ALLAA"
                # print mutations_string
                dummy = dummy + 1

                template = self.get_template(mutations_string)

                self.write_to_file( template, dummy )

            elif( len(self.protein_index) > 0.0):

                last_mutations_string = last_mutations_string + self.protein_index.pop(0)+" "+self.chain+" ALLAA\n\t"
                last = True


        if(last == True):        
            dummy = dummy + 1
            template = self.get_template(last_mutations_string)
            self.write_to_file( template, dummy )



        print "Length of protein is: ",self.protein_length

if __name__ == "__main__":
    run = InputFileFoldX()
    run.main()

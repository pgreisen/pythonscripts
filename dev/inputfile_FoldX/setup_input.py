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
                        chain = line[21:22]
                        resname = self.get_single_residue( line[17:20])

                        pdbstring = resname+chain+str(line[22:26]).strip()

                        self.protein_index.append( pdbstring )

        self.protein_length = len(self.protein_index)

        #self.min = min(self.protein_index)
        #self.max = max(self.protein_index)

    def get_template(self,mutations_string):

        template = '''<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>'''+self.pdbfile+''';
<BATCH>list.txt;
<COMMANDS>FOLDX_commandfile;
<PositionScan>#,'''+mutations_string+'''
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>false;
<pdb_hydrogens>false;
<complex_with_DNA> false;
<END>#;
<JOBEND>#;
<ENDFILE>#;'''
        return template



    def write_to_file( self, template, dummy ):
        with open(str(dummy)+"_run.txt",'w') as f:
            for line in template:
                f.write( line)



    def main(self):

        parser = argparse.ArgumentParser(description=" Generates input file for FoldX to compute mutational scan ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="pdbfile", help="Initial starting pdb")

        parser.add_argument("-p", "--processors", dest="processors", help="The number of processors to split the jobs into")

        input_variables = parser.parse_args()

        self.pdbfile = input_variables.pdbfile

        self.get_pdb_start_end()

        dummy = 0

        split_into = self.protein_length / int(input_variables.processors)

        print "Length of protein is: ",self.protein_length

        for i in range(1, self.protein_length+5, split_into):

            if( (len( self.protein_index )/split_into) > 0.0 ):

                mutations_string = self.protein_index.pop(0)+"a,"+self.protein_index.pop(0)+"a,"+self.protein_index.pop(0)+"a,"+self.protein_index.pop(0)+"a,"+self.protein_index.pop(0)+";"
                dummy = dummy + 1

                template = self.get_template(mutations_string)

                self.write_to_file( template, dummy )

            else:

                print self.protein_index, self.protein_length

                mutations_string = ""

                for i in self.protein_index:

                    mutations_string = mutations_string+","+str(i)+"a,"

                mutations_string[1:]



        print "Length of protein is: ",self.protein_length

if __name__ == "__main__":
    run = InputFileFoldX()
    run.main()

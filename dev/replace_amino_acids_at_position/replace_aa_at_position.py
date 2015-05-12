import sys
import sys, os, argparse

class ReplaceAA:

    def __init__(self):
        """
        @return:
        """
        self.position = -1
        self.aa = ""
        self.INSERT = True
        self.new_pdb = []
        self.input_pdb = ""
        self.offset = 0

    def insert_aa(self):

        with open(self.input_pdb,'r') as f:
            for line in f:


                if( str(line[22:26]).strip() == self.position and line[13:16] == "N" ):
                    new_line = line[0:17]+self.aa+line[20:]
                    self.new_pdb.append(new_line)

                elif( str(line[22:26]).strip() == self.position and str(line[13:16]).strip() == "CA" ):
                    new_line = line[0:17]+self.aa+line[20:]
                    self.new_pdb.append(new_line)

                elif( str(line[22:26]).strip() == self.position and str(line[13:16]).strip() == "C" ):
                    new_line = line[0:17]+self.aa+line[20:]
                    self.new_pdb.append(new_line)

                elif( str(line[22:26]).strip() == self.position and str(line[13:16]).strip() == "O" ):
                    new_line = line[0:17]+self.aa+line[20:]
                    self.new_pdb.append(new_line)


                elif( str(line[22:26]).strip() == self.position):
                    continue

                else:
                    self.new_pdb.append( line )



    def write_pdbfile(self):
        with open("sub_aa_"+self.aa+"_"+str(int(self.position)-self.offset)+".pdb", 'w') as f:
            for line in self.new_pdb:
                f.write(line)


    def main(self):

        parser = argparse.ArgumentParser(description="Replace an amino acid at a certain position - only inserted are backbone atoms (N CA C O)")
        # get the initial rosetta design as input
        parser.add_argument("--position", dest="position", help="Data set 1" )

        parser.add_argument("--aa", dest="aa", help="Insert this type of amino acid" )

        parser.add_argument("-f", "--pdbfile", dest="input_pdb", help="Input pdbfile" )

        input_variables = parser.parse_args()

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])


        self.insert_aa()
        self.write_pdbfile()


if __name__ == "__main__":
    run = ReplaceAA()
    run.main()

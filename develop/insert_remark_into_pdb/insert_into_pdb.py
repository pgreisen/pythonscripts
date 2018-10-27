import argparse, subprocess,os, shutil

# The script is executed in the following way:
# -l : three letter name of ligand residue
# -a : 3 letter code for amino acid 
# -r : residue number
# -n : cst number
# -f : file
# after inserting the new string into the file it will move
# the temporary file into the old file.
# The script is executed in the following manner:
#
#  python insert_into_pdb.py -f 2RF1_chainA_0001.pdb -l VXM -a TYR -r 100 -n 1
#  python insert_into_pdb.py -f 2RF1_chainA_0001.pdb -l VXM -a TYR -r 100 -n 2
#


class InsertRemarkIntoPDBFile:


    def __init__(self):
        pass

    def insert_remark_line(self,filename,ligandname, aminoacid, residuenumber, cstnumber):
        remarkline = "REMARK 666 MATCH TEMPLATE X "+str(ligandname)+"    0 MATCH MOTIF A "+str(aminoacid)+"  "+str(residuenumber)+"  "+str(cstnumber)+"  1\n"

        filelist = []
        with open(filename,'r') as f:
            for line in f:
                filelist.append( line )

        filelist.insert( int(cstnumber)-1,remarkline)
        
        tmpfile = open("tmp.txt",'w')
        for line in filelist:
            tmpfile.write(line)
        tmpfile.close()
        shutil.move("tmp.txt", filename)
        print "Remark line has been inserted into file"

    def main(self):

        parser = argparse.ArgumentParser(description="Insert remark into pdb file")

        parser.add_argument("-f", "--file", dest="inputfile",
                      help="File where remark will be inserted")

        parser.add_argument("-l", "--ligandname", dest="ligandname",
                      help="Name of ligand to be inserted in remark string")

        parser.add_argument("-a", "--amino_acid", dest="aa",
                      help="Amino acid for the cst")

        parser.add_argument("-r", "--residue_number", dest="nr",
                      help="Residue number of the amino acid")

        parser.add_argument("-n", "--cst_number", dest="cst_nr",
                      help="Cst number and which position it will be inserted into the pdb file")

        input_variable = parser.parse_args()

        self.insert_remark_line(input_variable.inputfile,input_variable.ligandname, input_variable.aa, input_variable.nr, input_variable.cst_nr)


if __name__ == "__main__":
   run = InsertRemarkIntoPDBFile()
   run.main()

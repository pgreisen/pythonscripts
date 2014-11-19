import argparse, subprocess,os, shutil

######################################################
#
#
# Combines files into
#
# python combinations.py
#
# the cst files are set in the main() function
#
######################################################


class CombineCstFiles:

    def __init__(self):
        self.position_one = ["LYS.cst","ARG.cst"]
        self.position_two = ["PP_2_CD1.cst","PP_2_CE1.cst","PP_2_CG.cst","TRP_2.cst","PP_2_CD2.cst","PP_2_CE2.cst","PP_2_CZ.cst","tstack_f_2_new.cst"]
        self.position_three = ["PP_1_CD1.cst","PP_1_CE1.cst","PP_1_CG.cst","TRP_1.cst","PP_1_CD2.cst","PP_1_CE2.cst","PP_1_CZ.cst","tstack_f_1_new.cst"]
        self.position_four = ["HB_O16.cst"]
        self.position_five = ["HB_O15.cst"]
        #self.position_six = ["HB_O10.cst"]

    def get_files(self,filename):
        with open(filename,'r') as f:
            tmp_file = f.readlines()
            tmp_file.append("\n")
        return tmp_file


    def write_cst_file(self,p1,p1_name,p2,p2_name,p3,p3_name,p4,p4_name,p5,p5_name,p6,p6_name):
        tmp_tot = p1+p2+p3+p4+p5+p6
        filename = p1_name+'_'+p2_name+'_'+p3_name+'_'+p4_name+'_'+p5_name+'_'+p6_name+".cst"

        with open(filename, 'w') as f:
            for line in tmp_tot:
                f.write(line)

    def combine_per_position(self):

        for i in self.position_one:
            p1_file = self.get_files( i )

            for j in self.position_two:
                p2_file = self.get_files( j )

                for k in self.position_three:
                    p3_file = self.get_files( k )

                    for l in self.position_four:
                        p4_file = self.get_files( l )

                        for m in self.position_five:
                            p5_file = self.get_files( m )

                            self.write_cst_file(p1_file,i.split('.')[0],p2_file,j.split('.')[0],p3_file,k.split('.')[0],p4_file,l.split('.')[0],p5_file,m.split('.')[0] )




    def main(self):

        self.combine_per_position()


    def test(self):
        a= combinations("1234")
        for i in a:
            print i


if __name__ == "__main__":
    run = CombineCstFiles()
    run.main()

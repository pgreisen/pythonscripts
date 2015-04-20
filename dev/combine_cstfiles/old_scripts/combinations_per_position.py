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
        self.position_one = ["PP_3_CD1.cst","PP_3_CD2.cst","PP_3_CE1.cst","PP_3_CE2.cst","PP_3_CG.cst","PP_3_CZ.cst","TRP_3.cst","tstack_f_3_new.cst"]
        self.position_two = ["PP_4_CD1.cst","PP_4_CD2.cst","PP_4_CE1.cst","PP_4_CE2.cst","PP_4_CG.cst","PP_4_CZ.cst","TRP_4.cst","tstack_f_4_new.cst" ]
        self.position_three = ["PP_5_CD1.cst","PP_5_CD2.cst","PP_5_CE1.cst","PP_5_CE2.cst","PP_5_CG.cst","PP_5_CZ.cst","TRP_5.cst","tstack_f_5_new.cst" ]



    def get_files(self,filename):
        with open(filename,'r') as f:
            tmp_file = f.readlines()
        return tmp_file


    def write_cst_file(self,p1,p1_name,p2,p2_name,p3,p3_name,):
        tmp_tot = p1+p2+p3
        filename = p1_name+'_'+p2_name+'_'+p3_name+".cst"

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
                    self.write_cst_file(p1_file,i.split('.')[0],p2_file,j.split('.')[0],p3_file,k.split('.')[0])


    def main(self):

        self.combine_per_position()


    def test(self):
        a= combinations("1234")
        for i in a:
            print i


if __name__ == "__main__":
    run = CombineCstFiles()
    run.main()
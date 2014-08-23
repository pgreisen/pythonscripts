import argparse, subprocess,os, shutil

######################################################
#
#
# Combines files into
#
# python combinations.py -c w90.cst,w205.cst,w43.cst
#
#
#
######################################################


class CombineCstFiles:

    def __init__(self):
        self.cstcollector = {}


    def combine_cst_files(self, list_of_cst_files):
        yield ''
        # (0, thing[0]), (1, thing[1]), (2, thing[2]), and so forth
        for i, d in enumerate( list_of_cst_files ):
            # being called recursively
            for comb in self.combine_cst_files( list_of_cst_files[i+1:]):
                yield d +"_" +comb


    def add_cst_files_to_list(self, cstfile ):
        coordinates = []
        with open(cstfile,'r') as f:
            for line in f:
                coordinates.append(line)
        coordinates.append("\n")
        self.cstcollector[cstfile] = coordinates


    def write_cst_to_files(self,cstfilename ):
        outputname = cstfilename[:-1]
        with open(outputname,'w') as f:
            tmpcstfile = cstfilename.split('_')
            for line in tmpcstfile:
                if( len(line) > 1):
                    tmpcst = self.cstcollector[line]
                    for newline in tmpcst:
                        f.write(newline)
        print "files is written",cstfilename



    def main(self):


        parser = argparse.ArgumentParser(description="Make combination of input cst-files")

        parser.add_argument("-c", "--cstfiles", dest="cstfiles",
                      help="cst files seperated by comma")

        input_variable = parser.parse_args()

        cstfiles = input_variable.cstfiles.split(',')

        # all possible combinations are made from the input cst-files
        tmp = self.combine_cst_files( cstfiles )

        # Loop over all the pairs
        for i in tmp:
            # list with cst files
            tmpcst = i.split('_')

            cstfilename = ""
            # skip any empty places in the list
            if( len( tmpcst ) > 2):
                for j in tmpcst:
                    if(len(j) > 2):
                        cstfilename += j+"_"
                        self.add_cst_files_to_list(j)
                self.write_cst_to_files( cstfilename )


    def test(self):
        a= combinations("1234")
        for i in a:
            print i


if __name__ == "__main__":
    run = CombineCstFiles()
    run.main()
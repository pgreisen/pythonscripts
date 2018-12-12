import sys, shutil, os, subprocess, argparse,re

'''
Argument1 is number of dockings to setup
Argument2 is path to input pdd

'''


class RemoveCBCB:


    def __init__(self):
        self.datafile = ""
        self.list_w_cutoffs = []

    def set_positions(self):

        with open(self.datafile,'r') as f:
            for line in f:
                tmpline = line.split(',')
                pairs = (tmpline[0],tmpline[1])
                self.list_w_cutoffs.append( pairs )


    def remove_files(self):
        dst = "Removed"
        try:
            shutil.os.mkdir( dst )
        except:
            print "Directory Removed already exists"

        files = os.listdir('.')

        for i in files:
            if( os.path.isfile( i ) and i.endswith(".fasta") ):

                for tmp in self.list_w_cutoffs:

                    tmp1 = tmp[0].split('_')[0]
                    tmp1_1 = tmp[0].split('_')[1]

                    tmp2 = tmp[1].split('_')[0]
                    tmp2_1 = tmp[1].split('_')[1]

                    re_str1 = tmp1+"._"+tmp1_1
                    re_str2 = tmp2+"._"+tmp2_1

                    search_group1 = re.search( r'^'+re_str1, i, re.M|re.I)
                    search_group2 = re.search( r'.*'+tmp2+'.?'+tmp2_1 , i , re.M|re.I)

                    if(search_group1 and search_group2):
                        shutil.copy(i , dst)
                        os.remove(i)
                        break

                    search_group1 = re.search( r'.*'+re_str1,  i, re.M|re.I)
                    search_group2 = re.search( r'^'+re_str2, i, re.M|re.I)

                    if(search_group1 and search_group2):
                        shutil.copy(i , dst)
                        os.remove(i)
                        break

                    search_group1 = re.search( r'.*'+re_str1,  i, re.M|re.I)
                    search_group2 = re.search( r'.*'+re_str2, i, re.M|re.I)

                    if(search_group1 and search_group2):
                        shutil.copy(i , dst)
                        os.remove(i)
                        break


    # Parameters to run rosetta on hyak with patchdock output
    # locate all parameter files in one central director
    def main(self):

        parser = argparse.ArgumentParser(description="Setup jobs protein-protein docking on Davinci - location of PDB file has to be full path ")

        parser.add_argument('-f',dest='datafile', help='Datafile with positions to remove')

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_positions()
        self.remove_files()




if __name__ == "__main__":
    run = RemoveCBCB()
    run.main()

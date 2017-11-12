import sys, shutil, os, subprocess, argparse

'''

Requires a directory with pdb files. Loops over all the pdbs and makes a directory
and submit the job.

Requires: User needs to set the path to the different parameters such as parameters, flags, xml, database, and executables. 

and is simply executed like this:

python submit_glide.py

'''


class GetSCResults:


    def __init__(self):
        self.score_sc = {}
        self.score_delta_ife = {}
        self.exe = {}
        self.init_name = "rmsd_ife.csv"
        self.exe = "/work/greisen/Projects/LigandBinding/VitaminD/Appendix/4_docking_designs/secondRoundTHCNonProNonDis_ForPer/run.sh "
        self.parameterfile = "/work/greisen/Projects/LigandBinding/VitaminD/14072015_THC_docking_random/EnzdesSetDockingFromNBR/THC/THC.fa.params"
        self.scoreterm = "SC"

    def set_sc(self,pdbfile):
        with open(pdbfile+"_0001.pdb",'r') as f:
            for line in f:
                if(len(line) < 2):
                    continue
                elif( line[0:2] == "AT"):
                    continue
                elif( line[0:2] == "HE"):
                    continue
                elif( line[0:2] == self.scoreterm ):
                    tmp = line.split()
                    self.score_sc[pdbfile] = tmp[1]
                else:
                    continue


    def write_to_file(self,directory_name):
        with open("dir_"+directory_name+"_ife_sc.dat",'w') as f:
            for key in self.score_delta_ife:
                f.write(key+","+self.score_sc[key]+","+self.score_delta_ife[key]+"\n")
        print "DONE"

    def main(self):
        parser = argparse.ArgumentParser(description="Rescore top 20 to get shape complementarity and dump a dat-file with numbers and names")
        # get the initial rosetta design as input
        parser.add_argument("--init_name", dest="init_name", help="Name of file to identify scores after first round (Defaul=rmsd_ife.csv )",default="rmsd_ife.csv" )
        input_variables = parser.parse_args()
        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])
        PTH = os.path.abspath('./')
        # Get the files in the directory
        files = os.listdir(PTH)
        for fl in files:
            if( os.path.isdir(fl) ):
                os.chdir( fl )
                # reset dictionary
                self.score_sc = {}
                self.score_delta_ife = {}

                # get file and assert it exists
                name_of_file = str(fl)+"_"+self.init_name
                if( os.path.exists( name_of_file ) ):
                    with open(name_of_file,'r') as f:
                        for line in f:
                            tmp = line.split(',')
                            self.score_delta_ife[ tmp[0] ] = tmp[2].strip()
                    # setup and run ?

                    for key in self.score_delta_ife:
                        exe = self.exe+" "+key+".pdb"+" "+self.parameterfile
                        subprocess.Popen(exe,shell=True).wait()
                        self.set_sc(key)
                        print "Directory is "+fl+" and the pdbfile is "+key
                    # write to file
                    self.write_to_file(fl)


                else:
                    print "file is not present: ", fl

                os.chdir('../')


if __name__ == "__main__":
    run = GetSCResults()
    run.main()

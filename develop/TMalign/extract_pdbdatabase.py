import sys, shutil, os, subprocess,argparse
from collections import OrderedDict, defaultdict


class TMalignAnalysis:


    def __init__(self):
        self.parent_pdb = ""
        self.data = defaultdict(list)  #{}
        self.aligned_pdbs = OrderedDict()
        self.datafile ="TM.sup_all_atm_lig"
        self.tmscore_threshold = 0.50


    def analyse_data(self):
        tofile = 0
        present = 0
        final_ter = 0
        key = ""
        donotadd = "1"
        tmpkeys = []
        go = 0
        with open(self.datafile) as f:
            for line in f:
                tmpline = line.split()
                if(tmpline[0] == "REMARK"):
                    if( tmpline[1] == "Chain" ):

                        if( len(tmpline) == 4 ):
                            size = tmpline[3]
                        else:
                            size = tmpline[4]
                        tmp_key = tmpline[2].split(':')[1]
                        newkey = tmp_key.split('.')[0]
                        if(key in self.aligned_pdbs.keys() ):
                            donotadd = "0"
                        else:
                            # continue
                            self.aligned_pdbs[newkey] = []

                        tmpkeys.append(newkey)

                        if( tmpline[2].split(':')[0] == '1'):
                            datakey = newkey

                    if(  tmpline[1] == "Aligned" ):
                        aligned_length = tmpline[3].split(',')[0]
                        rmsd = tmpline[5].split(',')[0]
                        tmp_tmscore =  tmpline[6].split('=')[1].split(',')[0]
                        self.data[datakey].append( size )
                        self.data[datakey].append( tmp_tmscore )
                        self.data[datakey].append( aligned_length )
                        self.data[datakey].append( rmsd )
                        if( float(tmp_tmscore) >= self.tmscore_threshold ):
                            self.aligned_pdbs[newkey] = []
                            go = 1

                elif(tmpline[0] == "TER" and go != 0):
                    tofile = 1
                    final_ter = 1

                elif(tmpline[0] == "ATOM" and tofile == 0 and donotadd == "1" and go != 0 ):
                    self.aligned_pdbs[ tmpkeys[0] ].append(line)

                elif(tmpline[0] == "ATOM" and tofile == 1 and donotadd == "1" and go != 0):
                    self.aligned_pdbs[tmpkeys[1]].append(line)
        return go

    def write_to_file(self):
        for key in self.aligned_pdbs:
            if( len( self.aligned_pdbs[key] ) == 0 ):
                continue

            with open(key+"_aligned.pdb",'w') as f:
                # print "The key is: ",key
                for line in self.aligned_pdbs[key]:
                    f.write(line)



    def main(self):
        parser = argparse.ArgumentParser(description="Insert mutations into fasta file ")
        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-m", "--mutations", dest="mutation", help="New amino acid entity")

        parser.add_argument("-n", "--native", dest="native", help="Old amino acid entity")

        parser.add_argument("-p", "--position", dest="position", help="Position to mutate", type=int )

        parser.add_argument("-c", "--chain", dest="chain", help="Light or heavy chain" )

        PTH = os.path.abspath('./')

        self.parent_pdb = sys.argv[1]
        # Get the files in the directory
        #pdbfiles = os.listdir(PTH)
        PTH = "/z/insulin/users/pjug/databases/pdb/"
        pdbfiles = os.listdir("/z/insulin/users/pjug/databases/pdb/")

        for pdbfile in pdbfiles:
            if(os.path.isfile(PTH+pdbfile) and str(pdbfile[-3:]) == 'pdb'):
                # print pdbfile
                shutil.copy(PTH+pdbfile,'./')
                exe = "~/Programs/tmalign/TMalign "+pdbfile+" "+self.parent_pdb+" -o TM.sup"
                subprocess.Popen(exe,shell=True).wait()
                os.remove(pdbfile)
                go = self.analyse_data()
                if( go == 1 ):
                    # print go, go== 1
                    self.write_to_file()


        with open("dataanalysis.dat", 'w') as f:
            f.write("PDB ID,Length Query,TMscore,Align,RMSD\n")
            for key in self.data:
                f.write(key+',')
                for value in self.data[key]:
                    f.write(value+',')
                f.write("\n")



if __name__ == "__main__":
   run = TMalignAnalysis()
   run.main()

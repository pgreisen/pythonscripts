import sys, shutil, os, subprocess, argparse,re

'''
Argument1 is number of dockings to setup
Argument2 is path to input pdd

'''


class RemoveCBCB:


    def __init__(self):
        self.datafile = ""
        self.fastafile = ""
        self.list_w_cutoffs = []
        self.designs_before = 0
        self.designs_after = 0

    def set_positions(self):

        with open(self.datafile,'r') as f:
            for line in f:
                tmpline = line.split(',')
                pairs = (tmpline[0],tmpline[1])
                self.list_w_cutoffs.append( pairs )



    def get_fasta_files_and_filter_distance(self,fastafile):
        designs = {}
        with open(fastafile, 'r') as f:
            for line in f:
                if(line[0] == '>'):
                    key=line[1:].strip()
                    designs[key] = ""
                else:
                    designs[key] += line.strip()
        # my_dict.pop('key', None)
        # S102_A,L85_A,
        # >A55E_D18N_R77T_0
        self.designs_before = len(designs.keys() )

        for i in self.list_w_cutoffs:
            pair1 = i[0].split('_')[0]
            pair2 = i[1].split('_')[0]
            for key in designs.keys():
                if(key.find(pair1) != -1 and key.find(pair2) != -1):
                    print(key,pair1,pair2)
                    designs.pop(key, None)
        self.designs_after = len(designs.keys() )
        with open("designs_filtered_based_on_distance.fasta",'w') as f:
            for line in designs.keys():
                f.write('>'+key+'\n')
                f.write(designs[key]+'\n')
        print("Designs before removal: ",self.designs_before)
        print("Designs after removal: ",self.designs_after)
                    
    # Parameters to run rosetta on hyak with patchdock output
    # locate all parameter files in one central director
    def main(self):

        parser = argparse.ArgumentParser(description="Setup jobs protein-protein docking on Davinci - location of PDB file has to be full path ")

        parser.add_argument('-d',dest='datafile', help='Datafile with positions to remove')
        parser.add_argument('-f',dest='fastafile', help='Datafile with positions to remove')


        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_positions()
        self.get_fasta_files_and_filter_distance(self.fastafile)
        #self.remove_files()

if __name__ == "__main__":
    run = RemoveCBCB()
    run.main()

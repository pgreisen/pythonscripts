import os,sys, argparse
from numpy import *
from numpy import linalg as LA
from collections import defaultdict
from collections import defaultdict

class GeneratePositionFile:

    def __init__(self):
        self.one_chain_only = False
        self.pdbfile = ""
        self.pdb = []
        self.length_pdb = 0
        self.chains = ([])
        self.chain_positions = defaultdict(list)


    def set_pdbfile(self):
        chain = []

        previous_chain = 'dummy'

        with open(self.pdbfile,'r') as f:
            chains_length = 1

            for line in f:

                if( line[0:4] == "ATOM"):

                    chainid = line[21:22]
                    chain.append( chainid )

                    if( chainid != previous_chain ):
                        self.chain_positions[previous_chain] = chains_length
                        previous_chain = chainid

                    # count c-alphas
                    if(line[13:15] == "CA"):
                        chains_length = chains_length + 1

        # assure that the last chain is added to the hash table as well
        # if( chainid != previous_chain ):
        self.chain_positions[previous_chain] = chains_length

        self.chains = set(chain)

    def main(self):

        parser = argparse.ArgumentParser(description="Generate a position file to be used in connection with matching")
        # get the initial rosetta design as input
        parser.add_argument("-d","--distance", dest="distance", help="Maximum distance to be included in the search for Cbeta-Cbeta distance (Default 5 AA)", type=float)

        parser.add_argument("-b", "--bundles", dest="helical_bundle", help="Four chains helical bundle with four chains is set to true", action="store_true", default=False )

        parser.add_argument("-f", "--pdbfile", dest="pdbfile", help="The name of the pdb file", default=None, type=str )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_pdbfile()


        r = sorted(self.chains, key=lambda item: (int(item.partition(' ')[0])
                                   if item[0].isdigit() else float('inf'), item))
        pos = 1
        # get the position files:
        number = 0
        for key in r:
            print key,self.chain_positions[key]

            number = self.chain_positions[key]
            print number

            while( pos < number):
                print pos
                pos = pos + 1
                if( pos+5 <= number ):
                    with open("pos_"+str(pos)+".pos",'w') as f:
                        tmp_line = "N_CST 4\n1: "+str(pos)+"\n2: "+str(pos+3)+"\n3: "+str(pos+4)+"\n4: ALL\n"
                        f.write(tmp_line)

if __name__ == "__main__":
    run = GeneratePositionFile()
    run.main()
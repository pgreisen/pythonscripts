#!/usr/bin/python
from optparse import OptionParser
import argparse



class GetChainIds():


    def __init__(self):
        self.debug = False
        self.chainid = 'A'


    def get_chain_id(self,pdbfile):
        chain_ids = []
        with open(pdbfile, 'r') as f:
            for line in f:
                if( line[0:4] == "ATOM"): # or line[0:4] == "HETA" ):
                    chain_ids.append( line[21:22] )
        return set( chain_ids )


    def replace_chain_id(self, pdbfile):
        new_pdbfile = open( "holo.pdb", 'w')
        with open(pdbfile, 'r') as f:
            for line in f:
                if( line[0:4] == "ATOM"):
                    new_line = line[0:21]+str(self.chainid)+line[22:]
                    new_pdbfile.write( new_line )
                else:
                    new_pdbfile.write(line)
        print "chain_id replaced"

    def main(self):
        parser = OptionParser()
        parser = argparse.ArgumentParser(description='get the all the chain ids present in pdb file')
        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="pdbfile", help="PDB file")
        inputs = parser.parse_args()
        chains_ids = self.replace_chain_id( inputs.pdbfile )
        # print chains_ids


if __name__ == '__main__':
    run = GetChainIds()
    run.main()

from collections import OrderedDict,defaultdict


class CleanPDB:


    def __init__(self):
        self.chains = ""
        self.pdbfile_chains = defaultdict(list)

        self.het_chains = defaultdict(list)

        self.hetatoms_in_chain = []

        self.ignore_hets = []

    def get_chains(self,pdbfile):
        with open(pdbfile,'r') as f:
            for line in f:

                if(line[0:4] == "HET "):
                    tmp = line.split()

                    a = tmp[1]+"_"+tmp[2]

                    self.hetatoms_in_chain.append(a)

                if(line[0:4] == "ATOM" ):

                    self.pdbfile_chains[line[21:22]].append(line)

                elif( line[0:4] == "HETA" ):
                    for tmp in self.hetatoms_in_chain:

                        tm = tmp.split('_')

                        if( line[17:20] == tm[0] and line[21:22] == tm[1] ):
                            key = line[17:20]+"_"+line[21:22]

                            self.het_chains[key].append(line)

        return self.pdbfile_chains, self.het_chains


'''
Assert chain and

Downloading PDB structure '1A46'...
[('TYS', 'I'), ('NA', 'H'), ('NA', 'H'), ('00K', 'H')]
Downloading PDB structure '1A5G'...
[('TYS', 'I'), ('NA', 'H'), ('NA', 'H'), ('00L', 'H')]
Downloading PDB structure '1N51'...
[('01B', 'B'), ('NH2', 'B'), ('MN', 'A2001'), ('MN', 'A2002'), ('MN', 'A2003')]
Downloading PDB structure '1QDQ'...
[('074', 'A')]

'''


'''
    def get_chains(self,pdbfile):
        with open(pdbfile,'r') as f:
            for line in f:

                if(line[0:4] == "HET "):
                    tmp = line.split()

                    a = (tmp[1],tmp[2])

                    self.hetatoms_in_chain.append(a)

                if(line[0:4] == "ATOM" or line[0:4] == "HETA"):

                    self.pdbfile_chains[line[21:22]].append(line)

        return self.pdbfile_chains, self.hetatoms_in_chain
'''
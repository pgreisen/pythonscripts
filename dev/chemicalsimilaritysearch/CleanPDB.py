from collections import OrderedDict,defaultdict


class CleanPDB:


    def __init__(self):
        self.chains = ""
        self.pdbfile_chains = defaultdict(list)



    def get_chains(self,pdbfile):
        with open(pdbfile,'r') as f:
            for line in f:
                if(line[0:4] == "ATOM" or line[0:4] == "HETA"):

                    import pdb; pdb.set_trace()

                    self.pdbfile_chains[line[21:22]].append(line)

        return self.pdbfile_chains
from pylab import *
import sys,os,  argparse
from collections import defaultdict

class FilterDesignsScore:

    def __init__(self):
        self.mutations_positions = defaultdict(list)
        self.ife = 0.0
        self.sc = 0.0
        self.sasa = 0.0
        self.unsatisfied = 0.0
        # self.scores_design = {}
        self.design_thresholds = {
            "interfE" : -10.79,
            "sc" : 0.69,
            "sasa" : 0.5
        }
        self.scores = ["ddg_complex","interfE","pstat","rwl","sc","unsat"]
        self.designs = []
        
    def evaluate_designs(self,pdbfile):
        design_scores = {}
        with open(pdbfile,'r') as f:
            for line in f:
                if(line[0:4] == "ATOM" or line[0:4] == "HETA"):
                    continue
                elif(len( line.split() ) == 2  ):
                    key = line.split()[0]
                    value = line.split()[1]
                    if(key in self.scores ):
                        design_scores[key] = float( value )
        return design_scores
    def write_designs(self):
        with open("designs_above_threshold.sh",'w') as f:
            f.write("mkdir DesignsAboveThreshold\n")
            f.write("for pdb in ")
            for pdb in self.designs:
                f.write(pdb+" ")
            f.write(";\ndo\ncp $pdb DesignsAboveThreshold/;\ndone")

    def main(self):
        path = "./"
        dirs = os.listdir( path )
        # loop over PDBs
        for file in dirs:
            if( os.path.isfile( file ) and file.endswith(".pdb") ):
                design_scores = self.evaluate_designs(file)
                sum = 0.0
                for key in design_scores:
                    if(key in self.design_thresholds):
                        if(design_scores[key] < 0 ):
                            if(design_scores[key] < self.design_thresholds[key]):
                                sum += 1
                        if(design_scores[key] > 0 ):
                            if(design_scores[key] > self.design_thresholds[key]):
                                sum += 1
                if(sum > 1.99):
                    self.designs.append(file)

        self.write_designs()





if __name__ == "__main__":
    run = FilterDesignsScore()
    run.main()

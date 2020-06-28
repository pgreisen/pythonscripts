from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import glob
import sys, getopt
from Bio import SeqIO, pairwise2
import Bio.SubsMat.MatrixInfo as matrices
import sklearn.cluster as cluster
import argparse
from Bio.SubsMat import MatrixInfo as matlist
import itertools
import multiprocessing
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from sklearn.neighbors.nearest_centroid import NearestCentroid
import random
random.seed(42)
'''
good example here:

https://stackoverflow.com/questions/42013822/parallelizing-four-nested-loops-in-python?rq=1

'''

class ClusterPDBs():

    def __init__(self):

        self.clusters = 15
        self.gap_open = -10
        self.gap_extend = -0.5
        self.seqs = {}
        self.seqs_list = []
        self.scores = [[]]
        self.fastafile = ""
        #self.scorematrixmatrix = matlist.blosum62
        self.scorematrixmatrix = matlist.pam30
        self.records = None
        

    def global_alignment(self,params):
        i = params[0]
        j = params[1]
        a = pairwise2.align.globalds(self.records[i].seq, self.records[j].seq, self.scorematrixmatrix, self.gap_open, self.gap_extend)
        # just writing it out for now
        (s1, s2, score, start, end) = a[0]
        return (s1, s2, score, start, end,i,j)


    def set_and_write_clusters(self):

        handle = open(self.fastafile, "rU")
        self.records = list(SeqIO.parse(handle, "fasta"))

        lookup = {}
        for i in self.records:
            lookup[i.seq] =  i.name

        nm_seqs = len(self.records)
        scores = [[0 for i in range(nm_seqs)] for j in range(nm_seqs)]

        # Generate values for each parameter
        seqA = range(nm_seqs)
        seqB = range(nm_seqs)

        # Generate a list of tuples where each tuple is a combination of parameters.
        # The list will contain all possible combinations of parameters.
        paramlist = list(itertools.product(seqA,seqB))
        ##paramlist = list(itertools.combinations(seqA,2))

        # Generate processes equal to the number of cores
        pool = multiprocessing.Pool()
        print("Number of cpus used: ",multiprocessing.cpu_count())
        # Distribute the parameter sets evenly across the cores
        # Global alignment returns this: (s1, s2, score, start, end,i,j)
        res = pool.map(self.global_alignment, paramlist)

        # dataframe to collect all the data
        columns = ["S1", "S2", "Score", "Start", "End", "i","j"]
        dfs = []
        for tmp in res:
            score_ = tmp[2]
            dict_ = {"S1" : tmp[0],
                     "S2" : tmp[1],
                     "Score" : score_,
                     "Start" : tmp[3],
                     "End" : tmp[4],
                     "i" : lookup[tmp[0]],
                     "j" : lookup[tmp[1]]}
            scores[tmp[-2]][tmp[-1]] = score_
            dfs.append(pd.DataFrame([dict_]))

        df = pd.concat(dfs)
        df = df.reset_index()
        df.to_csv("phylogenetic_data.csv")

        from sklearn.metrics import pairwise_distances_argmin_min
        from sklearn.cluster import AgglomerativeClustering 
        hc = AgglomerativeClustering(n_clusters = self.clusters, linkage='ward') 

        y_hc=hc.fit_predict(scores)
        cluster_labels = hc.labels_
               
        clf = NearestCentroid()
        clf.fit(scores, y_hc)

        target_counts = pd.Series(y_hc).value_counts()
        target_counts.plot.barh(colors=['#cb181d', '#fb6a4a', '#fcae91', '#fee5d9'],edgecolor='white');
        plt.title('Cluster Counts')
        plt.xlabel('Count')
        plt.ylabel('Cluser');
        plt.savefig('variants_per_clusters.png')        
        
        labels = hc.labels_
        clusters = [[] for i in range(self.clusters)]
        for i in range(0, nm_seqs):
            clusters[labels[i]].append(self.records[i])

        with open("cluster_centers_stochastic.fasta",'w') as f:
            for i in range(0, self.clusters):
                rd_ = random.choice(clusters[i])
                seq_ = rd_.seq
                name_ = rd_.name+"_cluster_"+str(i)
                f.write(">"+name_+"\n")
                f.write(str(seq_)+"\n")

        for i in range(0, len(clusters)):
            output_handle = open("c." + str(i) + ".fasta", "w")
            SeqIO.write(clusters[i], output_handle, "fasta")
            output_handle.close()

    def write_fasta_to_file(self):
        with open(self.tmpfastafile,'w') as f:
            for line in self.seqs:
                f.write(">"+line+'\n')
                f.write(str(self.seqs[line])+'\n')

    def main(self):
        parser = argparse.ArgumentParser(description="Cluster information")
        parser.add_argument('-g', dest='gap_extend', help='Gap extension', default=-0.5, type=float)
        parser.add_argument('-f', dest='fastafile', help='Fasta file with all sequence to be clustered')

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_and_write_clusters()


if __name__ == "__main__":
        run = ClusterPDBs()
        run.main()

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

        nm_seqs = len(self.records)
        scores = [[0 for i in range(nm_seqs)] for j in range(nm_seqs)]

        # Generate values for each parameter
        seqA = range(nm_seqs)
        seqB = range(nm_seqs)

        # Generate a list of tuples where each tuple is a combination of parameters.
        # The list will contain all possible combinations of parameters.
        paramlist = list(itertools.product(seqA,seqB))

        # Generate processes equal to the number of cores
        pool = multiprocessing.Pool()
        print("Number of cpus used: ",multiprocessing.cpu_count())
        # Distribute the parameter sets evenly across the cores
        res = pool.map(self.global_alignment, paramlist)

        for tmp in res:
            scores[tmp[-2]][tmp[-1]] = tmp[2]

        from sklearn.metrics import pairwise_distances_argmin_min
        kmeans = cluster.KMeans(self.clusters)

        results = kmeans.fit(scores)
        closest, _ = pairwise_distances_argmin_min(results.cluster_centers_, scores, metric="hamming")
        centroids = results.cluster_centers_

        labels = results.labels_
        clusters = [[] for i in range(self.clusters)]
        for i in range(0, nm_seqs):
            clusters[labels[i]].append(self.records[i])

        with open("cluster_centers.fasta",'w') as f:
            for i in range(len(closest)):
                seq_ = str(self.records[closest[i]].seq)
                name_ = self.records[closest[i]].name+"_cluster_"+str(i)
                f.write(">"+name_+"\n")
                f.write(seq_+"\n")

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

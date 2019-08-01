from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import glob
import sys, getopt
from Bio import SeqIO, pairwise2
import Bio.SubsMat.MatrixInfo as matrices
import sklearn.cluster as cluster
import argparse


class ClusterPDBs():

    def __init__(self):

        self.clusters = 9
        self.scorematrix = "blosum62"
        self.gap_open = -10
        self.gap_extend = -0.5
        self.seqs = {}
        self.seqs_list = []
        self.scores = [[]]
        self.fastafile = ""

    def set_and_write_clusters(self):

        handle = open(self.fastafile, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        from Bio.SubsMat import MatrixInfo as matlist
        matrix = matlist.blosum62
        nm_seqs = len(records)
        scores = [[0 for i in range(nm_seqs)] for j in range(nm_seqs)]

        for i in range(0, nm_seqs ):
            for j in range(0, nm_seqs ):
                # parallelize loop
                '''
                from multiprocessing import Pool
                pool = Pool(processes=4)  # start 4 worker processes
                result_list = pool.map(process_x, my_list)
                '''

                a = pairwise2.align.globalds(records[i].seq, records[j].seq, matrix, -10, -0.5)
                (s1, s2, score, start, end) = a[0]
                scores[i][j] = score

        kmeans = cluster.KMeans(self.clusters)
        results = kmeans.fit(scores)
        labels = results.labels_
        clusters = [[] for i in range(self.clusters)]
        for i in range(0, nm_seqs):
            clusters[labels[i]].append(records[i])

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

import boto3
import paramiko,time
import argparse
from boto3.s3.transfer import S3Transfer

class AddSignalPeptide:

    def __init__(self):
        self.protein_seqs = {}
        self.signal_peptide = ""
        self.fastafile = ""

    def add_signalpeptide_fasta_seqs(self):
        with open(self.fastafile,'r') as f:
            for line in f:
                if (line[0] == ">"):
                    key = line.strip()[1:]
                    self.protein_seqs[key] = self.signal_peptide
                else:
                    self.protein_seqs[key] += line.strip()
        self.write2file(self.protein_seqs)

    def write2file(self,dictionary):
        filename = self.fastafile.split(".")[0]
        with open("AddSignalPeptide_"+filename+".fasta", 'w') as f:

            for key in self.protein_seqs.keys():
                print(key)
                f.write(">"+key+"\n")
                f.write(self.protein_seqs[key]+"\n")


    def main(self):

        parser = argparse.ArgumentParser(description="Add a signal peptide")
        parser.add_argument('-s', dest='signal_peptide', help='Fasta file with signal peptide')
        parser.add_argument('-f', dest='fastafile', help='File containing the fasta sequences')

        args_dict = vars(parser.parse_args())

        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.add_signalpeptide_fasta_seqs()


if __name__ == "__main__":
    run = AddSignalPeptide()
    run.main()
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse, re
from Bio.Restriction import *
from Bio.Seq import Seq
from datetime import date
from Bio.Alphabet import generic_dna, generic_protein
# Class written to insert substitution
import InsertSubstitutions
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


class IDGFragments:


# prices taken from here (Tubes):
# https://www.idtdna.com/pages/products/genes-and-gene-fragments/gblocks-gene-fragments

    idt_gblock = {

        250: 71.10,
        500: 71.10,
        750: 80.10,
        1000: 89.10,
        1250: 143.10
    }


    def __init__(self):
        self.today = date.today()
        self.ref_dna_file = ""
        self.designs_aa_fasta_file = ""
        self.codon_table_file = ""
        self.cluster = 12
        self.cluster_size = 12
        self.ofile = str(self.today)+"_data"

        # parameters needed for the run
        # bsaI sites for cloning
        self.preseq = "GTCACGGTCTCA"
        self.postseq = "CGAGACCAGTCA"
        # range of basepair offset
        self.bp_offset_nterm = 3
        self.bp_offset_cterm = 3
        self.bp_offset = 3
        # variable parameter for the gene length
        # g-blocks to synthetize
        self.cutoff = [250, 500, 750, 1000]
        self.ref_dna = ""
        self.ref_key = ""




    def get_fasta_sequence( self,designs_aa_fasta_file):
        '''

        :param designs_aa_fasta_file: fasta file where header contains mutations in the format A55H
        :return: dictionary with header as key and sequence as value
        '''
        sequences_ = {}
        with open(designs_aa_fasta_file, 'r') as f:
            for line in f:
                if (line[0] == ">"):
                    tmp = line.replace('>','').split()
                    ids = tmp[0]
                    sequences_[ids] = ""
                else:
                    tmp = line.strip()
                    sequences_[ids] = sequences_[ids] + tmp
        return sequences_


    def insert_codon_mutation_into_DNA(self, sequences):
        '''

        :param sequences: contains correct named variants
        :return: mutated dictionary with DNA sequences
        '''

        # intialize method here:
        ris = InsertSubstitutions.InsertSubstitutions()
        # Insert substitutions into DNA string
        dna_seqs = {}
        for key in sequences.keys():
            stop_of_list = -1
            mut_header = ""
            tmpheader = key.split('_')
            integer_ = isinstance(tmpheader[-1], int)
            if (integer_ == False):
                stop_of_list = len(tmpheader)
            for i in tmpheader[0:stop_of_list]:
                mut_header += i + "_"
            fastaheader, sequence = ris.run(self.ref_dna_file, self.codon_table_file, mut_header[:-1])
            assert len(sequence) == len(self.ref_dna[self.ref_key])
            dna_seqs[fastaheader] = sequence
        return dna_seqs


    def get_price( self,length):
        ln = int(length)
        if ln <= 250:
            return 71.10
        elif ln <= 500:
            return 71.10
        elif ln <= 750:
            return 80.10
        elif ln <= 1000:
            return 89.10
        elif ln <= 1250:
            return 143.10
        elif ln <= 1500:
            return 179.1
        else:
            return "Error"


    def get_bsaI_analysis_and_write_fasta_file_w_variants(self, df):
        '''
        run bsaI analysis on DNA sequences and inserts results into dataframe
        :param df: dataframe with values
        :return: dataframe with bsaI analysis
        '''
        lengt = []
        with open(self.ofile+".fasta",'w') as f:
            for i, varname, dnaseq, st, end in zip(df.index, df["VarName"],df["seq(DNA)"],df["start_cut_DNA"],df["end_cut_DNA"]):
                f.write(varname+"\n")
                lengt.append(len(self.preseq+dnaseq[st:end]+self.postseq))
                idt_seq_ = self.preseq+dnaseq[st:end]+self.postseq
                a = Analysis(restrictionbatch=RestrictionBatch(['BsaI']),sequence=Seq(idt_seq_) )
                prot_sites_ = a.full()
                prot_key = list(prot_sites_.keys())[0]
                cutsites = ', '.join(map(str, prot_sites_[prot_key]))
                df.iloc[i, df.columns.get_loc("bsaI_sites")] = cutsites
                df.iloc[i, df.columns.get_loc("IDT_gene")] = idt_seq_
                f.write(idt_seq_+"\n")
        return df


    def reassign_missing_variants(self, df):
        '''
        :param df:
        :return:
        '''
        reassign_cluster = {}
        # sub-divide the rest into exciting cluster
        for i,j,k,c,l in zip(df.index, df["start_DNA_w_padding"], df["end_DNA_w_padding"], df["Clustered"],df["Size_of_gene"]):
            if(c == 0):
                for key in dict_w_start_stop:
                    if( (j >= dict_w_start_stop[key][0]) & (k <= dict_w_start_stop[key][1])):
                        # add and sort such that we can add it to the
                        # shortest gene cluster
                        tmpkey = "index_"+str(i)
                        reassign_cluster[tmpkey] = (i,j,k,c,key,dict_w_start_stop[key][2])
        # sort dictionary
        listofTuples = sorted(reassign_cluster.items() ,  key=lambda x: x[-1] )
        r_c = dict(listofTuples)
        for key in reassign_cluster.keys():
            df.iloc[r_c[key][0], df.columns.get_loc("Clustered")] = 1
            df.iloc[r_c[key][0], df.columns.get_loc("Cluster_name")] = r_c[key][4]
            df.iloc[r_c[key][0], df.columns.get_loc("Size_of_gene")] = r_c[key][5]
        return df


    def set_residue_numbers(self, seq_df):
        '''
        :param self:
        :param seq_df:
        :return:
        '''
        # initialize new columns
        seq_df["start_DNA"] = 0
        seq_df["end_DNA"] = 0
        # loop over keys in hash
        # Contain all the positions mutated
        residue_nr = []
        # contains residue numbers for each variant
        tmpresiduenr = {}
        for index, varname in zip(seq_df.index, seq_df["VarName"]):
            new_array = []
            nametmp_ = varname.split("_")
            for j in nametmp_:
                new_array.append( int(j[1:-1] ) )
            mn = min(new_array)
            mx = max(new_array)
            # The -3 is due to counting from 0 and protein vs DNA
            seq_df.iloc[index,seq_df.columns.get_loc("start_DNA")] = 3*mn-3
            #
            seq_df.iloc[index,seq_df.columns.get_loc("end_DNA")] = 3*mx
        return seq_df


    def set_clusters(self, df, init_values, cutoff):
        for i, j, l, c in zip(df.index, df["start_DNA_w_padding"], df["length_DNA_w_padding"], df["Clustered"]):

            # bool value so you dont overwrite existing clusters
            if (c == 1):
                continue
            elif ((j in init_values) & (l < cutoff)):
                df.iloc[i, df.columns.get_loc("Clustered")] = 1
                df.iloc[i, df.columns.get_loc("Cluster_name")] = "cluster_" + str(j) + "_" + str(j + cutoff) + "_" + str(
                    cutoff)
                df.iloc[i, df.columns.get_loc("Size_of_gene")] = cutoff
            else:
                continue
        return df


    def setup_dataframe(self, dna_seqs):
        '''

        :param self:
        :return:
        '''
        # Create a dataframe with the AA and DNA sequences
        df_sequences = pd.DataFrame.from_dict(dna_seqs, orient='index', columns=['seq(DNA)'])
        df_sequences.reset_index(inplace=True)
        df_sequences = df_sequences.rename(columns={'index': 'VarName'})
        df_sequences["VarName"] = df_sequences["VarName"].str.replace(">", "")

        df_sequences = self.set_residue_numbers(df_sequences)
        df_sequences["start_DNA_w_padding"] = df_sequences["start_DNA"] - self.bp_offset
        df_sequences["end_DNA_w_padding"] = df_sequences["end_DNA"] + self.bp_offset

        df_sequences.sort_values(by=["start_DNA_w_padding", "end_DNA_w_padding"], inplace=True)
        df_sequences["length_DNA_w_padding"] = df_sequences["end_DNA_w_padding"] - df_sequences["start_DNA_w_padding"]

        df_sequences["Clustered"] = 0
        df_sequences["Cluster_name"] = ""
        df_sequences["Size_of_gene"] = 0

        df = df_sequences.copy()
        df = df.reset_index()

        return df

    def cluster_seqs_and_write_to_file(self,df):
        '''
        :param self:
        :return:
        '''
        df_cluster = DataFrame(df, columns=['start_DNA_w_padding', 'end_DNA_w_padding'])
        kmeans = KMeans(n_clusters=self.cluster_size).fit(df_cluster)
        centroids = kmeans.cluster_centers_
        #plt.scatter(df_cluster['start_DNA_w_padding'], df_cluster['end_DNA_w_padding'], c=kmeans.labels_.astype(float), s=50,alpha=0.5)
        #plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)

        cluster_map = pd.DataFrame()
        cluster_map['data_index'] = df.index.values
        cluster_map['cluster'] = kmeans.labels_

        # merge dataframes
        df_kmeans = pd.merge(right=df, right_index=True, left=cluster_map, left_index=True)

        df_kmeans.drop(labels=['data_index', 'index', \
                               'start_DNA_w_padding', 'end_DNA_w_padding', 'length_DNA_w_padding', \
                               'Clustered', 'Size_of_gene', 'Cluster_name'], axis=1, inplace=True)
        nr_clusters = list(set(kmeans.labels_))
        df_kmeans['start_cut_DNA'] = 0
        df_kmeans['end_cut_DNA'] = 0
        df_kmeans['Size_of_gene'] = 0
        df_kmeans['Cluster_name'] = ""

        for cluster_ in nr_clusters:
            # get each cluster
            subcluster = df_kmeans[df_kmeans['cluster'] == cluster_]
            # start cut DNA
            mn_ = min(subcluster['start_DNA'])
            # end cut DNA
            mx_ = max(subcluster['end_DNA'])
            # size of gene
            length_ = mx_ - mn_
            # cluster name
            clustername_ = "Cluster_" + str(mn_) + "_" + str(mx_) + "_cluster_nr_" + str(cluster_)
            for i in df_kmeans.index:
                if (df_kmeans.iloc[i, df_kmeans.columns.get_loc('cluster')] == cluster_):
                    df_kmeans.iloc[i, df_kmeans.columns.get_loc('start_cut_DNA')] = mn_ - self.bp_offset
                    df_kmeans.iloc[i, df_kmeans.columns.get_loc('end_cut_DNA')] = mx_ + self.bp_offset
                    df_kmeans.iloc[i, df_kmeans.columns.get_loc('Size_of_gene')] = length_ + self.bp_offset
                    df_kmeans.iloc[i, df_kmeans.columns.get_loc('Cluster_name')] = clustername_

        # Calculate price
        df_kmeans["IDT_price"] = df_kmeans["Size_of_gene"].apply(self.get_price)
        print("The price of the order is: ", df_kmeans["IDT_price"].sum())
        print("Number of unique cluster: ", len(df_kmeans["Cluster_name"].unique()))
        tmp = df_kmeans.groupby("Cluster_name").count()
        print("Size of the different clusters: ", tmp["VarName"])
        # sets the IDT gene and bsaI cut-sites
        df_kmeans["bsaI_sites"] = ""
        df_kmeans["IDT_gene"] = ""
        df_kmeans = self.get_bsaI_analysis_and_write_fasta_file_w_variants( df_kmeans)

        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna, generic_protein
        df_kmeans["AA_seq"] = ""
        for i, j in zip(df_kmeans.index, df_kmeans["IDT_gene"]):
            tmpseq_ = j[len(self.preseq):-len(self.postseq)]
            my_seq = Seq(tmpseq_)
            df_kmeans.iloc[i, df_kmeans.columns.get_loc("AA_seq")] = my_seq.translate()
        df_kmeans.to_excel(self.ofile + ".xlsx", index=False)
        return df_kmeans


    def main(self):
        # pass argument to class
        parser = argparse.ArgumentParser(description="")
        # get the initial rosetta design as input
        parser.add_argument("-r", "--ref", dest="ref_dna_file", help="Reference DNA sequence to insert mutations into")
        parser.add_argument("-d", "--designed_amino_acids_file", dest="designs_aa_fasta_file",help="Contains the designs as strings e.g. >A34K_H88M")
        parser.add_argument("-c", "--codontable", dest="codon_table_file",help="Contains the codon table to be used for the insertions")
        parser.add_argument("-s", "--clustersize", dest="cluster", help="Number of clusters the sequences are divided into",default=12, type=int)



        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # get reference DNA
        self.ref_dna = self.get_fasta_sequence(self.ref_dna_file)
        # get the AA sequence of the designs
        sequences = self.get_fasta_sequence(self.designs_aa_fasta_file)
        self.ref_key = list(self.ref_dna.keys())[0]

        # Insert substitution into DNA sequence
        dna_sequences = self.insert_codon_mutation_into_DNA( sequences )


        # setup data into dataframe
        df = self.setup_dataframe(dna_sequences)

        # cluster sequences and return
        #df = self.set_clusters(df)
        df = self.cluster_seqs_and_write_to_file( df )



if __name__ == "__main__":
   run = IDGFragments()
   run.main()




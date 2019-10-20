#!/usr/bin/env python
import pandas as pd
import argparse
import datetime
import numpy as np

class VerifyNameSeq:
    """
    Take a fasta-file with multiple sequences and their header name e.g. >A12F_G36K with a sequence
    and verifies that the sequence and name corresponds.
    The script is executed

    python VerifyNameSeq.py -f FASTAFILE

    if there is a mismatch between name and sequence it return a report with the problematic names.

    """

    def __init__(self):
        self.fastaseq = {}
        self.fastafile = ""
        self.character_to_split_on = "_"
        self.dataframe_seqs = None
        self.name = ""
        self.today = str(datetime.date.today())

    def generate_dataframe(self):
        self.dataframe_seqs = pd.DataFrame.from_dict(self.fastaseq, orient='index').reset_index().rename(columns={'index': 'Name', 0: 'sequence'})

    def set_fastasequence_from_file(self ):
        '''
        :param fastafile
        :return: sets dictionary with names and sequences
        '''

        with open(self.fastafile, 'r') as f:
            for i in f:
                if (i[0] == ">"):
                    key = i.strip()[1:]
                    self.fastaseq[key] = ""
                else:
                    self.fastaseq[key] += i.strip()
        # generate a dataframe which can be used as report to the user
        self.generate_dataframe()

    def check_if_name_matches(self):
        '''

        :return: dictionary with mismatches
        '''
        report = {}
        for i in self.fastaseq.keys():
            muts = i.split(self.character_to_split_on)
            report_ = ""
            for j in muts:
                wt_ = j[0]
                # index in python starts at 0
                pos_ = int(j[1:-1])-1
                aa_ = j[-1]
                if( aa_ !=  self.fastaseq[i][pos_]):
                    report_ += wt_+str(j[1:-1])+self.fastaseq[i][pos_]+","
            report[i] = report_[:-1]
        return report


    def highlight_greaterthan(self, s, columnname="Mismatch"):

        if s[columnname] != '':
            return ['background-color: red']*3
        else:
            return ['background-color: white']*3


    def get_report(self,report):

        # merge report to input file
        df_report = pd.DataFrame.from_dict(report, orient='index').reset_index().rename(columns={'index': 'Name', 0: 'Mismatch'})
        dftot = pd.merge(right=self.dataframe_seqs,right_on="Name",left=df_report,left_on="Name" )

        # highlight_table = dftot.style.apply(self.highlight_last_row)
        # highlight_table = dftot.style.applymap(color_negative_red)
        highlight_table = dftot.style.apply(self.highlight_greaterthan, axis=1)

        df_html = highlight_table.render()
        with open(self.today+"_mismatch_table_highlighted.html",'w') as f:
            for line in df_html:
                f.write(line)

        dftot.to_html(self.today + "_" + self.name + ".html",index=False)
        dftot.to_html(self.today+"_"+self.name+".csv",index=False)

    def main(self):
        parser = argparse.ArgumentParser(description="Verifies if names and sequences match - e.g. A14F does it contain F at position 14")
        parser.add_argument("-f", "--file", dest="fastafile", help="Fasta file with names of the sequence")
        parser.add_argument("-c", "--character", dest="charac", help="Fasta file with names of the sequence",default="_")
        parser.add_argument("-n", "--name", dest="name", help="Name of output file with date added to it", default="test")


        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_fastasequence_from_file()
        report = self.check_if_name_matches()
        self.get_report(report)


if __name__ == "__main__":
   run = VerifyNameSeq()
   run.main()

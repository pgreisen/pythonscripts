import argparse, subprocess, os
import collections


class GetGermline:
    def __init__(self):
        self.query_sequence = ""
        # Score, E-value, Pct_align,Seq_cov
        self.germlines = collections.OrderedDict() # {}
        self.query = collections.OrderedDict() #{}
        self.germline_file = ""

        self.germlines_data = ["VH.dat","VJ.dat"]
        self.path = "./"

        self.VH = "/Users/pgreisen/pythonscripts/dev/Igblast_analysis/Human/20160428_human_VH.fasta"
        self.VJ = "/Users/pgreisen/pythonscripts/dev/Igblast_analysis/Human/20160428_human_VJ.fasta"

        self.germline_vj_ids = []
        self.germline_vh_ids = []
        self.germline_vd_ids = []
        self.top_hits = 3
        self.vh_seq = collections.OrderedDict() #{}
        self.vj_seq = collections.OrderedDict() #{}


    def set_query_sequence(self):
        with open(self.germline_file,'r') as f:
            for line in f:
                if(line[0] == '>'):
                    key = line[1:].strip()
                    self.query[key] = ""
                else:
                    self.query[key] += self.query[key]+line.strip()


    def get_germlines_reference(self,germlinefile):
        sequences = {}
        with open(germlinefile, 'r') as f:
            for line in f:
                if (line[0] == ">"):
                    # tmp = line.split()
                    ids = line[1:].strip() # tmp[0]
                    sequences[ids] = ""
                else:
                    tmp = line.strip()
                    sequences[ids] = sequences[ids] + tmp
        return sequences

    def get_ig_sequence(self, variable_part):
        if(variable_part == "VH.dat"):
            vh = self.get_germlines_reference(self.VH)
            for key in self.germlines.keys():
                for key2 in vh.keys():
                    if(key == key2):
                        self.vh_seq[key] = vh[key2]

        elif (variable_part == "VJ.dat"):
            vj_seq = {}
            vj = self.get_germlines_reference(self.VJ)
            for key in self.germlines.keys():
                for key2 in vj.keys():
                    if (key == key2):
                        self.vj_seq[key] = vj[key2]


    def set_germline_ids(self,igblastoutput):
        with open(igblastoutput,'r') as f:
            dummy = 0
            dummy2 = 0
            for line in f:
                if( len(line) > 3 and line[0:3] == "lcl"  and dummy <= self.top_hits):
                    tmpline = line.split()
                    key = tmpline[0].split('|')[1].strip()
                    if(key not in self.germlines.keys() ):
                        self.germlines[key] = ""
                    self.germlines[key] = self.germlines[key]+tmpline[-2]+","+tmpline[-1]+","
                    dummy += 1

                if(len(line) > 10 and line[0:10] == "Alignments"):
                    tmp_bool = True

                if ( line[0] == "V" and tmp_bool == True and dummy2 <= self.top_hits):
                    tmpline = line.split()
                    key = tmpline[3].strip()
                    self.germlines[key] = self.germlines[key]+tmpline[1]+","+tmpline[2]
                    dummy2 += 1

        self.get_ig_sequence(igblastoutput[-6:])

    def write_to_file(self,name,variable_domain):
        with open(name+'.fasta','w') as f:
            for v_key in variable_domain:
                f.write(">" + v_key + "\n")
                f.write(variable_domain[v_key] + "\n")

    def write_table_to_file(self):
        with open("IG_table.csv",'w') as f:
            f.write("Germline,Score,E-value,Pct_align,Seq_cov\n")
            for line in self.germlines.keys():
                f.write(line+','+self.germlines[line]+'\n')

    def main(self):
        parser = argparse.ArgumentParser(description='Get reference germline')
        parser.add_argument("-f", "--file", dest="germline_file", help="Germline file")
        parser.add_argument("-s", "--query", dest="query_sequence", help="The query sequence")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        dirs = os.listdir(self.path)
        for filename in dirs:
            if( filename[-6:] in self.germlines_data ):
                self.set_germline_ids( filename )

        self.write_to_file("VJ_germlines", self.vj_seq)
        self.write_to_file("VH_germlines", self.vh_seq)

        self.write_table_to_file()


if __name__ == "__main__":
   run = GetGermline()
   run.main()


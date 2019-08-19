import sys
import argparse, csv


class DiffFasta:

    def __init__(self):

        self.fastafile1 = ""
        self.fastafile2 = ""
        self.all_designs = ""
        self.fasta_1 = {}
        self.fasta_2 = {}
        self.offset = 0
        self.mutations = []
        self.key_a = ""
        self.key_b = ""
        self.design_seq = {}
        self.prefix = "CSA_"
        self.sciworm_format = {}

    def set_fasta_file(self, fastafile, container):
        # dummy variable in case fasta names are not unique
        dummy = 0
        with open(fastafile, 'r') as f:
            for line in f:
                if (line[0] == ">"):
                    key = line.strip()[1:] + str(dummy)
                    container[key] = ""
                    dummy += 1
                    print(dummy)
                else:
                    container[key] = container[key] + line.strip()

    def sequence_a_b(self, a, b):
        i = 1
        # print "A is ", a
        for key in a:
            self.key_a = key
            seq_a = a[key]
            seqlengthA = len(a[key])

        seq_b = b
        seqlengthB = len(b)

        assert seqlengthA == seqlengthB

        dummy = 1 + self.offset
        diff_dummy = ""
        for i in range(seqlengthA):

            if (seq_a[i] != seq_b[i]):
                diff_dummy = diff_dummy + seq_a[i] + str(dummy) + seq_b[i] + "_"
                self.mutations.append(seq_a[i] + str(dummy) + seq_b[i])
            dummy += 1
        return diff_dummy

    def set_designs(self, datafile, sequences):
        dummy = 1
        with open(datafile, 'r') as f:
            for line in f:
                if (line[0] == ">"):
                    # 2016-12-07 change such that it goes for group id
                    # tmp = line.split()
                    # ids = tmp[0]
                    tmp = line.strip().split("_")

                    # ids = tmp[-1]
                    ids = line[1:].strip() + "_" + str(dummy)
                    dummy += 1

                    if (ids not in sequences.keys()):
                        sequences[ids] = {}
                        sequences[ids]["LC"] = ""
                        sequences[ids]["HC"] = ""
                else:
                    tmp = line.strip()
                    # Todo fix this towards generic algorithm
                    if (tmp[0] == "M"):
                        sequences[ids]["LC"] = sequences[ids]["LC"] + tmp
                    else:
                        print("Debug code")

    def main(self):
        parser = argparse.ArgumentParser(description=" ")
        parser.add_argument("-s", "--file1", dest="fastafile1", help="Native light chain", default="")
        parser.add_argument("-p", "--file2", dest="fastafile2", help="Native heavy chain", default="")
        # parse file with sequences
        parser.add_argument("-a", "--all", dest="all_designs", help="Fasta format")
        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_designs(self.all_designs, self.design_seq)

        if (len(self.fastafile1) != 0):
            self.set_fasta_file(self.fastafile1, self.fasta_1)

        if (len(self.fastafile2) != 0):
            self.set_fasta_file(self.fastafile2, self.fasta_2)

        for key in self.design_seq.keys():
            mutations = ""

            if (len(self.design_seq[key]["HC"]) != 0 and len(self.design_seq[key]["LC"]) != 0):
                hc_mutations = self.sequence_a_b(self.fasta_2, self.design_seq[key]["HC"])
                mutations = hc_mutations + "HC_"
                lc_mutations = self.sequence_a_b(self.fasta_1, self.design_seq[key]['LC'])
                mutations = mutations + lc_mutations + "LC"
                new_key = self.prefix + mutations + "_" + key.split("_")[-1]
                self.sciworm_format[new_key] = (self.design_seq[key]["HC"], self.design_seq[key]["LC"])

            elif (len(self.design_seq[key]["HC"]) != 0):

                hc_mutations = self.sequence_a_b(self.fasta_2, self.design_seq[key]["HC"])
                mutations = hc_mutations + "HC"
                new_key = self.prefix + mutations  # +"_"+key.split("_")[-1]
                self.sciworm_format[new_key] = (self.design_seq[key]["HC"], "")

            elif (len(self.design_seq[key]["LC"]) != 0):

                lc_mutations = self.sequence_a_b(self.fasta_1, self.design_seq[key]['LC'])
                mutations = mutations + lc_mutations
                new_key = self.prefix + mutations  # +"_"+key.split("_")[-1]
                self.sciworm_format[new_key] = ("", self.design_seq[key]["LC"])
            else:
                print("Nothing to do!!!!!!")
                continue

        with open("sciworm.fasta", 'w') as f:
            for key in self.sciworm_format.keys():
                if (self.sciworm_format[key][0] != ""):
                    f.write(">" + key + "\n")
                    f.write(self.sciworm_format[key][0] + "\n")
                if (self.sciworm_format[key][1] != ""):
                    f.write(">" + key + "\n")
                    f.write(self.sciworm_format[key][1] + "\n")


if __name__ == "__main__":
    run = DiffFasta()
    run.main()

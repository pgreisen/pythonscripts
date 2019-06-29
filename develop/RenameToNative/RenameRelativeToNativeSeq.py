import sys
import argparse,csv

class DiffFasta:
    '''
    import sys
    sys.path.append(PATH)
    import RenameRelativeToNativeSeq as rrn
    r_ = rrn.DiffFasta()
    var_ = r_.diff_sequence_a_b(nn_seq_,test_seq_)

    '''



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
        self.prefix = "0365-0001-1767_"
        self.postfix = ""
        self.chain_lc = ""
        self.chain_hc = ""
        self.sciworm_format = {}
        self.keep = 0

    def set_fasta_file(self,fastafile,container):

        with open(fastafile,'r') as f:
            for line in f:
                if( line[0] == ">"):
                    key = line.strip()[1:]
                    container[key] = ""
                else:
                    container[key] = container[key]+line.strip()
        # print container,key

    def sequence_a_b(self, a, b):
        '''
        :param a: dictionary with keys and value as sequence
        :param b: template sequence
        :return: mutations between them
        '''
        i = 1
        # print "A is ", a
        for key in a:
            print(key)
            self.key_a = key
            seq_a = a[key]
            seqlengthA = len(a[key])

        seq_b = b
        seqlengthB = len(b)
        print("B is: ", b)
        # print key
        print(seqlengthA, seqlengthB, a, b)

        assert seqlengthA == seqlengthB

        dummy = 1 + self.offset
        diff_dummy = ""
        for i in range(seqlengthA):

            if (seq_a[i] != seq_b[i]):
                # print seq_a[i], seq_b[i]
                # print key_a,key_b
                diff_dummy = diff_dummy + seq_a[i] + str(dummy) + seq_b[i] + "_"
                self.mutations.append(seq_a[i] + str(dummy) + seq_b[i])
            dummy += 1
        return diff_dummy

    def diff_sequence_a_b(self, seq_a, seq_b):
        '''
        :param a: fasta sequence
        :param b: template sequence
        :return: mutations between them
        '''
        i = 1
        seqlengthA = len(seq_a)
        seqlengthB = len(seq_b)
        ##print "B is: ", seq_b
        # print key
        ##print seqlengthA, seqlengthB, seq_a, seq_b

        assert seqlengthA == seqlengthB

        dummy = 1 + self.offset
        diff_dummy = ""
        for i in range(seqlengthA):
            if (seq_a[i] != seq_b[i]):
                diff_dummy = diff_dummy + seq_a[i] + str(dummy) + seq_b[i] + "_"
                self.mutations.append(seq_a[i] + str(dummy) + seq_b[i])
            dummy += 1
        return diff_dummy

    def set_designs(self,datafile,sequences):

    	with open(datafile,'r') as f:
    	    for line in f:
                if(line[0] == ">"):
                    # 2016-12-07 change such that it goes for group id
                    tmp = line.strip().split("_")
                    ids = line[1:].strip()
                    if(ids not in sequences.keys()):
                        sequences[ids] = {}
                        sequences[ids]["LC"] = ""
                        sequences[ids]["HC"] = ""
                else:
                    tmp = line.strip()
                    if (self.chain_lc == "LC"):
                        sequences[ids]["LC"] = sequences[ids]["LC"] + tmp
                    elif (self.chain_hc == "HC"):
                        sequences[ids]["HC"] = sequences[ids]["HC"] + tmp
                    else:
                        print("Debug code or no chain was set for this run")

    def main(self):
        parser = argparse.ArgumentParser(description=" ")
        parser.add_argument("-s", "--file1", dest="fastafile1", help="Native light chain", default="")
        parser.add_argument("-p", "--file2", dest="fastafile2", help="Native heavy chain", default="")
        # PG 2017-03-22
        parser.add_argument("--prefix", dest="prefix", help="Native heavy chain", default="")
        parser.add_argument("--postfix", dest="postfix", help="ID to pair the two chains", default="")

        # parse file with sequences
        parser.add_argument("-a", "--all", dest="all_designs", help="Fasta format")

        # LC/HC
        parser.add_argument("--hc", dest="chain_hc", help="Fasta format")
        parser.add_argument("--lc", dest="chain_lc", help="Fasta format")

        parser.add_argument("--keep", dest="keep", help="Keep header from fasta file as part of name (Default: No)")

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_designs(self.all_designs,self.design_seq)

        if (len( self.fastafile1 ) != 0):
            self.set_fasta_file(self.fastafile1, self.fasta_1)

        if (len(self.fastafile2) != 0):
            self.set_fasta_file(self.fastafile2, self.fasta_2)

        for key in self.design_seq.keys():
            mutations = ""
            # add key to prefix
            if(self.keep != 0):
                self.prefix = key+"_"

            if(len(self.design_seq[key]["HC"]) != 0 and len(self.design_seq[key]["LC"]) != 0 ):
                hc_mutations = self.sequence_a_b(self.fasta_2, self.design_seq[key]["HC"])
                mutations = hc_mutations+"HC_"
                lc_mutations = self.sequence_a_b(self.fasta_1,self.design_seq[key]['LC'])
                mutations = mutations+lc_mutations+"LC"+"_"+str(self.postfix)
                new_key = self.prefix+mutations+"_"+key.split("_")[-1]+"_"+self.postfix
                self.sciworm_format[new_key] = (self.design_seq[key]["HC"],self.design_seq[key]["LC"])

            elif(len(self.design_seq[key]["HC"]) != 0):
                hc_mutations = self.sequence_a_b(self.fasta_2, self.design_seq[key]["HC"])
                mutations = hc_mutations+"HC"+"_"+str(self.postfix)
                new_key = self.prefix+mutations # +"_"+key.split("_")[-1]+"_"+self.postfix
                self.sciworm_format[new_key] = (self.design_seq[key]["HC"],"")

            elif( len(self.design_seq[key]["LC"]) != 0):
                lc_mutations = self.sequence_a_b(self.fasta_1,self.design_seq[key]['LC'])
                mutations = mutations+lc_mutations+"LC"+"_"+str(self.postfix)
                new_key = self.prefix+mutations # +"_"+key.split("_")[-1]+"_"+self.postfix
                self.sciworm_format[new_key] = ("",self.design_seq[key]["LC"])
            else:
                print("Nothing to do!!!!!!")
                continue

        with open("sciworm.fasta",'w') as f:
            for key in self.sciworm_format.keys():
                if( self.sciworm_format[key][0] != ""):
                    f.write(">"+key+"\n" )
                    f.write(self.sciworm_format[key][0]+"\n" )
                if( self.sciworm_format[key][1] != ""):
                    f.write(">"+key +"\n" )
                    f.write(self.sciworm_format[key][1]+"\n" )

if __name__ == "__main__":
   run = DiffFasta()
   run.main()

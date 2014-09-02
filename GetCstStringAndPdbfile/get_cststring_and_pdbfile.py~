from optparse import OptionParser
import shutil,os


def get_remark_string(file,aminoacids):
    with open(file,'r') as f:
        for line in f:
            if(line[0:4] == 'REMA'):
                for j in aminoacids:
                    if j in line:
                        substring = "1  1"
                        if substring in line:
                            shutil.copy(file,"substringmatch/")

def main():
    parser = OptionParser()

    parser.add_option('-f',dest='pdbfile',help='Pdbfile with remark line')
    parser.add_option('-a',dest='aminoacids',default=[],help='aminoacids',type=str)
    # note to store_true / store_false would not change the value. 
    parser.add_option('-c',dest='copytodir',action="store_true", help='Copy the files to directory - substringmatch',default=False)

    (options,args) = parser.parse_args()

    aminoacids = options.aminoacids.split(',')

    if(options.copytodir):

        # print os.path.exists("substringmatch")
        assert (os.path.exists("substringmatch")) == True

    get_remark_string(options.pdbfile,aminoacids)    




if __name__ == '__main__':
    main()

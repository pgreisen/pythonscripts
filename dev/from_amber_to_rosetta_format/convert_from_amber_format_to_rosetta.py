import sys
from pdbfile import *

'''

missing add accupancy to one instead of zero which is the output from amber

'''


class ConvertFromAmberFormatToRosetta:

    def main(self):
        file = sys.argv[1]
        deliminator = sys.argv[2]

        pd = pdbfile()
        # get pdbfile
        pdblist = pd.read_file(file)

        pd.split_into_multiple_pdb_files( pdblist, deliminator )


if __name__ == "__main__":
   run = ConvertFromAmberFormatToRosetta()
   run.main()
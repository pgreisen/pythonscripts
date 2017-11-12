from getligandposition import *

import argparse

class convert_to_alanines:


    def __init__(self):
        global pdbname
        pass


    """

    @requires:

    """

    def insert_alanine_at_position(self, pdbfile, positions):
        global pdbname

        pwrite = open(pdbname+"_alanine.pdb",'w')

        for line in pdbfile:

            if(line[0:4] == "ATOM" and line[23:26].strip() in positions ):

                if(line[12:15].strip() in ['N','CA','C','O','CB']):

                    new_line = line[0:17]+"ALA"+line[20:]
                    pwrite.write(new_line)
            else:
                pwrite.write(line)

        pwrite.close()

    # @requires: pdbfile, position file, and ligand name
    def main(self):

        parser = argparse.ArgumentParser(description='Narrows down position')
        parser.add_argument("-f", "--file", dest="pdbfile",
                      help="File with aligned ligand")
        input_var = parser.parse_args()

        gp = GetLigandPositions()

        # Get pdbfile
        p_pdbfile = gp.get_pdbfile(input_var.pdbfile)
        global pdbname
        pdbname = input_var.pdbfile

        # Get ligand coordinates
        ligand_coor = gp.get_ligandatom_coordinates( p_pdbfile )

        # return protein positions
        new_positions = gp.get_positions( p_pdbfile, ligand_coor)

        # debug
        print "The number of positions interacting with the ligand is ", len(new_positions), new_positions

        # Change amino acids at these positions to alanine #
        self.insert_alanine_at_position(p_pdbfile, new_positions)


if __name__ == "__main__":
   pos = convert_to_alanines()
   pos.main()
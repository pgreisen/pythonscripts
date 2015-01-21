import os,shutil,sys
from translate_rotate import * 
from pdbfile import * 
import os,sys, argparse

'''


'''


class AlignResidues:


    def __init__(self):
        self.targetresidues = {}
        self.queryresidues = {}
        self.targetresiduescoordinates = []
        self.queryresiduescoordinates = []

    # Returns coordinate array
    def set_array_from_dictionary(self):

        for key in self.targetresidues:

            if key in self.queryresidues :
                self.targetresiduescoordinates.append( self.targetresidues[ key ] )
                self.queryresiduescoordinates.append( self.queryresidues[ key ])

    def main(self):

        pdbfile_instance = pdbfile()




        parser = argparse.ArgumentParser(description="Solving a SVD for a number of residues")
        # get the initial rosetta design as input

        parser.add_argument("-f", "--targetfile", dest="targetfile", help="The target file for the superpositioning of the coordinates", default=None )


        parser.add_argument("-q", "--queryfile", dest="queryfile", help="The query file for the superpositioning of the coordinates", default=None )


        parser.add_argument("--targetresidues", dest="targetresidues", help="The residues used for the superpositioning", default=None )


        parser.add_argument("--queryresidues", dest="queryresidues", help="The query-residues used for the superpositioning", default=None )


        input_variables = parser.parse_args()


        # The pdb objects as lists
        targetfile = pdbfile_instance.read_file(input_variables.targetfile)
        queryfile = pdbfile_instance.read_file(input_variables.queryfile)


        # residues used for the superpositioning


        targetresidues = input_variables.targetresidues.split(',')
        queryresidues = input_variables.queryresidues.split(',')


        # def get_all_residueatoms(self, pdblist, residuename, residuenumber):


        for i in targetresidues:
            tmp = i.split()
            pdbfile_instance.get_all_residueatoms( targetfile, tmp[0], tmp[1], self.targetresidues )



        for j in queryresidues:
            tmp = j.split()
            pdbfile_instance.get_all_residueatoms( queryfile, tmp[0], tmp[1], self.queryresidues )


        # debug
        #print self.queryresidues
        #print self.targetresidues
        #print self.queryresidues
        #assert 1 == 2

        # Data transformation
        # Order matters hence list is given as argument
        self.set_array_from_dictionary()

        # Make translation / rotation
        tr_instance = translate_rotate()

        # debug 1
        tr_matrix, rot_matrix = tr_instance.get_rotate_translate(self.queryresiduescoordinates, self.targetresiduescoordinates )

        # debug 1
        ## print tr_matrix, rot_matrix



        # Apply translation / rotation
        query_all_atoms = pdbfile_instance.get_pdbcoordinates( queryfile )



        query_coor = tr_instance.transform_ligand_coordinates( query_all_atoms ,tr_matrix,rot_matrix)
        # print query_coor

        t_coordinates = pdbfile_instance.get_transformed_coordinates_pdbfile( query_coor, queryfile )


        ## print t_coordinates

        #####protein = pdbfile_instance.get_protein_coordinates(pose)


        filename = 'realign'+input_variables.queryfile
        pdbfile_instance.write_file( t_coordinates, filename )



if __name__ == "__main__":
    run = AlignResidues()
    run.main()

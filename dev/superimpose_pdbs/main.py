import os,shutil,sys
from translate_rotate import * 
from pdbfile import * 

'''


The script is executed like this 

python main.py pdb1.pdb pdb2.pdb

where it aligns pdb2 onto pdb1 using CA atoms. It is necessary that the two 
chains have the same length.

'''

class SuperPositionPDBs:

    def __init__(self):
        pass

    # Returns coordinate array
    def get_array_from_dictionary(self, list_of_atoms, dictionary):
        coordinate_array = []
        for i in list_of_atoms:
            coordinate_array.append(dictionary[i])
        return coordinate_array


    def main(self):

        pdbfile_instance = pdbfile()
        # Input
        # Write it as parser
        fl = sys.argv[1]

        # Get list of pdbfile
        query_pose = pdbfile_instance.read_file(fl)
        # Do the alignment per atom ID ( index of the pdb file )
        superposition_atom = "CA"
        query_atoms_coordinates = pdbfile_instance.get_atoms_per_name( query_pose, superposition_atom )

        # target pdb
        fl2 = sys.argv[2]
        target_pose = pdbfile_instance.read_file(fl2)
        # Input: listinstance, atoms_to_get, residue_of_these_atoms, residuename
        target_atoms_coordinates = pdbfile_instance.get_atoms_per_name(target_pose, superposition_atom )

        # debug
        # print target_atoms_coordinates
        # print len( target_atoms_coordinates ), len( query_atoms_coordinates )
        assert len( target_atoms_coordinates ) == len( query_atoms_coordinates )

        # Make translation / rotation
        tr_instance = translate_rotate()
        # superimpose query on target
        tr_matrix, rot_matrix = tr_instance.get_rotate_translate( query_atoms_coordinates, target_atoms_coordinates )

        # debug print tr_matrix
        # print tr_matrix
        # fine - there is a translation

        # getting the coordinates of all the atoms from the query_pose
        pdb_all_atoms = pdbfile_instance.get_all_pdbatoms( query_pose )

        # perform the transformation
        pdb_coor = tr_instance.transform_ligand_coordinates_lists(pdb_all_atoms, tr_matrix, rot_matrix)

        # debug print list with the array of each of the new coordinate pairs
        # print pdb_coor

        # getting the new coordinates from the transformation
        t_coordinates = pdbfile_instance.get_t_coordinates(pdb_coor, query_pose )

        protein = pdbfile_instance.get_protein_coordinates( query_pose )

        # Only aligned ligand
        all = t_coordinates

        filename = 'realign'+fl
        pdbfile_instance.write_file( all, filename)


if __name__ == "__main__":
    run = SuperPositionPDBs()
    run.main()

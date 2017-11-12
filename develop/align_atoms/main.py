import os,shutil,sys
from translate_rotate import * 
from pdbfile import * 

'''


The script is executed like this 

python main.py MPAT_1a53_nohet_1_relax_patchdock.out.1.pdb MPA_0001.fa.pdb

where it aligns MPA_0001.fa.pdb onto MPAT_1a53_nohet_1_relax_patchdock.out.1.pdb.
The parameters needed to be set are

# Atom names
query_atoms = ['P', 'O1','O2']
# Residue name
query_residuename = 'LIG'

# Atom names
target_atoms = ['P1', 'O2', 'O1']
# Residue name
target_residuename = '<0>'


'''


# Returns coordinate array
def get_array_from_dictionary(list_of_atoms,dictionary):

    coordinate_array = []
    for i in list_of_atoms:
        coordinate_array.append(dictionary[i])
        #    for value in dictionary.values():
        #        coordinate_array.append(value)
    return coordinate_array


def main():

    pdbfile_instance = pdbfile()

    # Input 
    # Write it as parser
    fl = sys.argv[1]

    # Get list of pdbfile
    pose = pdbfile_instance.read_file(fl)

    # Query
    # Atom names
    query_atoms = ['FE','NA','NB','NC','ND','C3B'] #,'O1']
    # Residue name
    query_residuename = 'HEM'

    # Input: listinstance, atoms_to_get, residue_of_these_atoms, residuename

    # Do the alignment per atom names
    query_atoms_coordinates = pdbfile_instance.get_atoms(pose,query_atoms,query_residuename)

    # Do the alignment per atom ID ( index of the pdb file )
    ##query_atoms = ['28','21','10'] #,'O1']
    ##query_atoms_coordinates = pdbfile_instance.get_atoms_per_index(pose,query_atoms,query_residuename)

    fl2 = sys.argv[2]
    target_pose = pdbfile_instance.read_file(fl2)
    # Target
    # Atom names
    target_atoms =  ['FE','NA','NB','NC','ND','C3B']
    # Residue name
    target_residuename = 'HEM'
    # Input: listinstance, atoms_to_get, residue_of_these_atoms, residuename
    # Return: dictionary
    target_atoms_coordinates = pdbfile_instance.get_atoms(target_pose,target_atoms,target_residuename)
    ##target_atoms_coordinates = pdbfile_instance.get_atoms_per_index(target_pose,target_atoms,target_residuename)

    # Data transformation
    # Order matters hence list is given as argument 
    query_coordinates = get_array_from_dictionary(query_atoms,query_atoms_coordinates)
    target_coordinates = get_array_from_dictionary(target_atoms,target_atoms_coordinates)

    # Debug the coordinates are different
    # print target_atoms_coordinates
    # print query_coordinates


    # Make translation / rotation
    tr_instance = translate_rotate()
    tr_matrix, rot_matrix = tr_instance.get_rotate_translate(target_coordinates,query_coordinates)


    # Apply translation / rotation
    ligand_all_atoms = pdbfile_instance.get_all_residueatoms(target_pose,target_residuename)
    
    ligand_coor = tr_instance.transform_ligand_coordinates(ligand_all_atoms,tr_matrix,rot_matrix)

    ####

    t_coordinates = pdbfile_instance.get_transformed_coordinates(ligand_coor,target_pose)

    protein = pdbfile_instance.get_protein_coordinates(pose)
    
    # Aligned system
    ##all = protein + t_coordinates

    # Only aligned ligand
    all = t_coordinates


    filename = 'realign'+fl2

    pdbfile_instance.write_file(all,filename)


if __name__ == "__main__":
    main()

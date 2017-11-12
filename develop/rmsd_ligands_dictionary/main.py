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

# @returns list with none hydrogen atoms
def get_ligand_atoms_none_hydrogen(pdbfile):
    ligand_atom_name = []
    ligand_name = ""
    for line in pdbfile:
        if(line[0:3] == "HET" and line[13] != 'H'):

            atom_name = line.split()[2]
            ligand_atom_name.append(atom_name)
            ligand_name = line.split()[3]
    return ligand_atom_name,ligand_name
    


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
    # NativePose
    fl = sys.argv[1]
    # DesignPose
    fl2 = sys.argv[2]

    
    # Get list of pdbfile
    native_pose = pdbfile_instance.read_file(fl)
    native_atoms,native_residuename = get_ligand_atoms_none_hydrogen(native_pose)

    # Input: listinstance, atoms_to_get, residue_of_these_atoms, residuename    
    native_atoms_coordinates = pdbfile_instance.get_atoms(native_pose,native_atoms,native_residuename)

    design_pose = pdbfile_instance.read_file(fl2)
    design_atoms, design_residuename = get_ligand_atoms_none_hydrogen(design_pose)
    
    # we need to sort the atoms because the order is different in some of them
    native_atoms.sort()
    design_atoms.sort()

    # Necessary
    # For debug
    # print native_atoms,design_atoms
    # outcomment 07-08-2015 for NAD+
    ##assert(native_atoms == design_atoms)
    # Return: dictionary
    design_atoms_coordinates = pdbfile_instance.get_atoms(design_pose,design_atoms,design_residuename)
    # Data transformation
    # Order matters hence list is given as argument 
    design_coordinates = get_array_from_dictionary(design_atoms,design_atoms_coordinates)
    native_coordinates = get_array_from_dictionary(native_atoms,native_atoms_coordinates)
    # Make translation / rotation
    tr_instance = translate_rotate()
    rmsd = tr_instance.get_dictionary_rmsd(native_atoms_coordinates, design_atoms_coordinates)
    ##rmsd = tr_instance.get_rmsd(design_coordinates,native_coordinates,design_atoms,native_atoms)

    print "The rmsd is ",fl2,round(rmsd,3)

if __name__ == "__main__":
    main()

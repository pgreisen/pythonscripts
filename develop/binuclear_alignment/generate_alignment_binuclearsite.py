#!/usr/bin/env python
import os,shutil
from get_binuclear_metalsite import *
from get_atoms_binuclear import *
from set_alignment import *
from write_binuclearpdb_files import *

# Requires PDB file
# Returns tmp.pdb, superimposed.pdb, sample_substrate.pdb, cry.pdb

# Initialization of classes
no = get_binuclear_metalsite()
ga = get_atoms_binuclear()
sa = set_alignment()
wpf = write_binuclearpdb_files()

# Alignments
HYDROXY = 0
HYDROXY_4atoms = 1
PHOPHORYL = 0
PHOSPHORYL_4atoms = 0
DISTANCE = 2.7

## Should be optimized
def write_aligned_coordinates(list_of_coordinates,name):
    wr = open('chi_'+str(name[0])+'.pdb','w')
    lngth = len(list_of_coordinates)
    while lngth > 0:
        crd = list_of_coordinates.pop(0)
        lngth = lngth - 1
        for i in crd:
            wr.write(i)
            
def get_vector_pdb(line):
    vec = Vector(float(line[31:39]),float(line[39:47]),float(line[47:55]))
    return vec

def generate_bridge(bridge_coordinates):
    template_string = 'HETATM 3008  O   HOH A 880     -27.684  28.326  44.634  1.00 17.62           O'
    x = wpf.set_length_digit(str(wpf.set_number_digits(bridge_coordinates[0])))
    y = wpf.set_length_digit(str(wpf.set_number_digits(bridge_coordinates[1])))
    z = wpf.set_length_digit(str(wpf.set_number_digits(bridge_coordinates[2])))
    bridge = template_string[0:30]+x+y+z+template_string[55:]
    return bridge

def generate_fraction(diff_vector,fraction):
    x = diff_vector[0]*fraction
    y = diff_vector[1]*fraction
    z = diff_vector[2]*fraction
    return Vector(x,y,z)


### MOdify here if empty return something!!!! write file bridgeatom.pdb
# @Requires list of hetero atoms,metal ions
# Returns bridge atom
def get_bridge_atom(het_list,metal_ion_1,metal_ion_2,protein_ligands):
    bridge = ''
    for j in het_list:
        tmp = get_vector_pdb(j)
        if (tmp - metal_ion_1).norm() < DISTANCE:
            if (tmp - metal_ion_2).norm() < DISTANCE:
                bridge = j
    if bridge == '':
        tmp = (metal_ion_1 + metal_ion_2)/2
        for i in protein_ligands:
            temp = generate_fraction((tmp - protein_ligands[i]),0.1)
            tmp = tmp + temp
        bridge = generate_bridge(tmp)
    return bridge

def write_bridge(bridge):
    output = open('bridgeatom.pdb','w')
    output.write(bridge)
    output.close()
    
# Write cry.pdb file
# @requires dictionary with metal ions, bridge_atom
# write metal ions, bridging oxygen, along with other het atoms
# if specified
def write_cry_monohydroxy_pdb(metal_ions,bridge_atom):
    crystal_coor = 'cry.pdb'
    fl = open(crystal_coor,'w')
    # change 05-08-2010 first originally 0 and second 1 to rotate
    fl.write(metal_ions[1])
    fl.write(metal_ions[0])
    fl.write(bridge_atom)
    fl.close()


def write_cry_pdb(metal_ions,bridge_atom,het_atm):
    crystal_coor = 'cry.pdb'
    fl = open(crystal_coor,'w')
    fl.write(metal_ions[0])
    fl.write(metal_ions[1])
    fl.write(str(bridge_atom).rstrip()+'\n')
    fl.write(het_atm)
    fl.close()

# @Requires substrate, which alignment to perform,
# PDB-file to get name, and make name of directory
def write_aligned_pdb_move_to_dir(substrate,PDB,VIZ,ID=0):
    atomic_coor = []
    crystal_coor = 'cry.pdb'
    test = str(PDB).split('.')
    path,templatefile = sa.set_path_substrate(substrate)
    # substrate_atoms, ligand name from a parameter file
    # Return parameter for evaluation of rmsd measure???????
    atomic_coor.append(sa.get_aligned_coor(path,VIZ,templatefile,crystal_coor))
    if atomic_coor[0] != 'RMSD TOO HIGH':
        write_aligned_coordinates(atomic_coor,test)
        # Making directory for rosetta input
        dir_name = test[0]+'_'+str(ID)
        os.mkdir(dir_name)
        # PDB dir_name
        shutil.move('chi_'+PDB,dir_name+'/chi_'+dir_name+'.pdb')
        shutil.copy('bridgeatom.pdb',dir_name+'/bridge_atom.pdb')
        # Move pdb file to directory with all pdb files
        print 'Done',PDB
        shutil.copy(PDB,dir_name+'/'+dir_name+'.pdb')
    else:
        print 'RMSD alignment too bad'


        
# Generate the alignement
def generate_alignment_all(VIZ,METAL,PDB):
    NO_BRIDGE = 1
    # Crystal coordinates
    crystal_coor = 'cry.pdb'
    # Getting list of transformed coordinates - append to a list 
    # PDB file to be aligned to - here make loop over all files
    file_object = no.get_obj_pdbfile(PDB)
    # ZN coordinates and heteroatoms are retrived
    metal_coor,het = no.get_metal_ions(file_object,METAL)
    # Get heteroatoms around metal ion
    a = no.get_hetero_around_metal(metal_coor,het)
    nn = {}
    # 'a' contains the metal ligands 
    p_ligg,b = no.get_protein_ligand_metal(file_object,metal_coor,nn)
    het_lst =[]
    if len(a) != 0:
        het_atom = a.keys()
        het_lst = ga.get_hetatoms_pdb(PDB,het_atom)
    # Operating with list 
    metal_ions = ga.get_metalion_pdb(PDB,METAL)
    # Getting the two metal ions
    metal_ion_1 = get_vector_pdb(metal_ions[0])
    metal_ion_2 = get_vector_pdb(metal_ions[1])
    # Returns list with hetatoms
    print 'het_atomhet_atomhet_atomhet_atom',het_atom
#    het_lst = ga.get_hetatoms_pdb(PDB,het_atom)
    bridge = get_bridge_atom(het_lst,metal_ion_1,metal_ion_2,p_ligg)
    write_bridge(bridge)
    # Alignment to choose
    # 0 hydroxy alignment, 1 Hydroxy_4 alignment
    # Dummy variable maybe it can be removed
    ID = 0
    if(HYDROXY and NO_BRIDGE == 1):
        substrate = 0
        ALIGNMENT = 0
        # Setting the alignment 
        # Taking from the naming of the ligand in its LG1_0001.pdb
        substrate_atoms = ['ZN1','ZN2','O4']
        # Write tmp.pdb file for alignment 
        sa.set_align_file(substrate_atoms,ALIGNMENT)
        write_cry_monohydroxy_pdb(metal_ions,bridge)
        write_aligned_pdb_move_to_dir(substrate,PDB,VIZ,ID)
    elif(HYDROXY_4atoms and NO_BRIDGE == 1):
        substrate = 1
        # Substrate atoms
        substrate_atoms = ['ZN1','ZN2','O4','O1']
        ALIGNMENT = 1
        sa.set_align_file(substrate_atoms,ALIGNMENT)
        # Removing bridge atom from list
        for elem in het_lst:
            if elem == bridge:
                het_lst.remove(bridge)
        number_of_crystal_hetatms = len(het_lst)
        for i in range(number_of_crystal_hetatms):
            write_cry_pdb(metal_ions,bridge,het_lst[i])
            ID = ID + 1
            write_aligned_pdb_move_to_dir(substrate,PDB,VIZ,ID)
    if(PHOPHORYL and NO_BRIDGE == 1):
        substrate = 2
        ALIGNMENT = 2
        # Setting the alignment 
        # Taking from the naming of the ligand in its LG1_0001.pdb
        substrate_atoms = ['ZN1','ZN2','O4']
        # Write tmp.pdb file for alignment 
        sa.set_align_file(substrate_atoms,ALIGNMENT)
        write_cry_monohydroxy_pdb(metal_ions,bridge)
        write_aligned_pdb_move_to_dir(substrate,PDB,VIZ,ID)
    elif(PHOSPHORYL_4atoms and NO_BRIDGE == 1):
        substrate = 3
        # Substrate atoms
        substrate_atoms = ['ZN1','ZN2','O4','O1']
        ALIGNMENT = 3
        sa.set_align_file(substrate_atoms,ALIGNMENT)
        for elem in het_lst:
            if elem == bridge:
                het_lst.remove(bridge)
        # print het_lst
        # Removing bridge atom from list
        # het_lst.remove(bridge)
        number_of_crystal_hetatms = len(het_lst)
        for i in range(number_of_crystal_hetatms):
            write_cry_pdb(metal_ions,bridge,het_lst[i])
            ID = ID + 1
            write_aligned_pdb_move_to_dir(substrate,PDB,VIZ,ID)

def main():
    # Constants for output file
    # For viz
    # VIZ = 'END\n'
    VIZ = 'TER\n'
    METAL = ['ZN','MG','MN','FE','CU','NI','CO','CA']
    #    path = './'
    PDB_file = sys.argv[1]
    # Getting the directory
    # origdir = os.getcwd()
    #    files = os.listdir(path)
    #    for PDB_file in files:
    #        if os.path.isfile(PDB_file):
    #    test = PDB_file.split('.')
    # If any files present without extension we test it here
    #    lngth = len(test)
    #    if lngth > 1 and test[1] == 'pdb':
    generate_alignment_all(VIZ,METAL,PDB_file)
                
                
if __name__ == "__main__":
    main()
    
    

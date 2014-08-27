#!/usr/bin/env python
import os,shutil
from numpy import *
from get_binuclear_metalsite import *
from get_atoms_binuclear import *
from set_alignment4atoms import *
from write_binuclearpdb_files import *

# Initialization of classes
no = get_binuclear_metalsite()
ga = get_atoms_binuclear()
sa = set_alignment()
wpf = write_binuclearpdb_files()


# Constant used to get ligand atoms and get bridge atom
DISTANCE = 2.7
# RMSD value for the alignment - empirical number
RMSD_THRESHOLD = 0.8
# Atoms from the ligand which should be used in the alignment
# The order is quite necessary to follow
SUBSTRATE_ATOMS = ['ZN1','ZN2','O5','O1']


# Should be optimized
def write_aligned_coordinates(list_of_coordinates,name):

    wr = open('chi_'+str(name[0])+'.pdb','w')
    lngth = len(list_of_coordinates)
    while lngth > 0:
        crd = list_of_coordinates.pop(0)
        lngth = lngth - 1
        for i in crd:
            wr.write(i)
            
def get_vector_pdb(line):
    vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
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
    return array([x,y,z])


### MOdify here if empty return something!!!! write file bridgeatom.pdb
# @Requires list of hetero atoms,metal ions
# Returns bridge atom
def get_bridge_atom(het_list,metal_ion_1,metal_ion_2,protein_ligands):
    bridge = ''

    for j in het_list:
        tmp = get_vector_pdb(j)
        if linalg.norm(tmp - metal_ion_1) < DISTANCE:
            if linalg.norm(tmp - metal_ion_2) < DISTANCE:
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
    

def write_cry_pdb(metal_ions,bridge_atom,het_atm):
    crystal_coor = 'cry.pdb'
    fl = open(crystal_coor,'w')
    fl.write(metal_ions[0])
    fl.write(metal_ions[1])
    # Hack 05-10-2010
    # fl.write(str(bridge_atom).rstrip()+'\n')
    # fl.write(het_atm)
    fl.write(het_atm)
    fl.write(str(bridge_atom).rstrip()+'\n')
    fl.close()

# @Requires substrate, which alignment to perform,
# PDB-file to get name, and make name of directory
def write_aligned_pdb_move_to_dir(LIGANDFILE,PDB,VIZ,ID=0):
    atomic_coor = []
    crystal_coor = 'cry.pdb'

    test = str(PDB).split('.')
        
    RMSD_BOOLEAN,atomic_coor = sa.get_aligned_coor(VIZ,LIGANDFILE,crystal_coor,RMSD_THRESHOLD,SUBSTRATE_ATOMS)


    if(RMSD_BOOLEAN):
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
# What goes on here?
def generate_alignment_all(VIZ,METAL,PDB,LIGANDFILE):

    # What does this do - magic constant?
    #    NO_BRIDGE = 1
    NO_BRIDGING_OXYGEN = 1

    
    # Crystal coordinates temporery file
    ## 16-08-2012 Commented out
    # crystal_coor = 'cry.pdb'
    
    # Getting list of transformed coordinates - append to a list 
    file_object = no.get_obj_pdbfile(PDB)

    # ZN coordinates and heteroatoms are retrived
    metal_coor,het = no.get_metal_ions(file_object,METAL)

    # Get heteroatoms around metal ion
    hetero_atoms_around_metal_ions = no.get_hetero_around_metal(metal_coor,het,DISTANCE)

    protein_ligands = no.get_protein_ligand_metal(file_object,metal_coor,DISTANCE) 

    list_of_hetero_atoms = []

    if len(hetero_atoms_around_metal_ions) != 0:
        het_atom = hetero_atoms_around_metal_ions.keys()
        list_of_hetero_atoms = ga.get_hetatoms_pdb(PDB,het_atom)

    # Operating with list 
    metal_ions = ga.get_metalion_pdb(PDB,METAL)

    # Getting the two metal ions
    metal_ion_1 = get_vector_pdb(metal_ions[0])
    metal_ion_2 = get_vector_pdb(metal_ions[1])
    
    # Returns list with hetatoms
    bridge = get_bridge_atom(list_of_hetero_atoms,metal_ion_1,metal_ion_2,protein_ligands)
    write_bridge(bridge)

    # Dummy variable
    ID = 0

    for elem in list_of_hetero_atoms:
        if elem == bridge:
            list_of_hetero_atoms.remove(bridge)

    number_of_crystal_hetatms = len(list_of_hetero_atoms)

    # 17-08-2012
    for i in range(number_of_crystal_hetatms):
        write_cry_pdb(metal_ions,bridge,list_of_hetero_atoms[i])
        ID = ID + 1
        write_aligned_pdb_move_to_dir(LIGANDFILE,PDB,VIZ,ID)


'''
Isolate main and make the program more readable

'''
def main():

    # Constants for output file
    # For viz
    # VIZ = 'END\n'

    VIZ = 'TER\n'
    METAL = ['ZN','MG','MN','FE','CU','NI','CO','CA']
    #    path = './'
    PDB_file = sys.argv[1]
    LIGAND_FILE = sys.argv[2]
    
    generate_alignment_all(VIZ,METAL,PDB_file,LIGAND_FILE)
                
                
if __name__ == "__main__":
    main()
    
    

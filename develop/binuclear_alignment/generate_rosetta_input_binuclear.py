#!/usr/bin/env python
'''
How to run the script and which inputs are needed

0. Main should take inputs and how to run program

1. generate_output_rosetta is impossible to read - make it set up the run

'''


from numpy import *
from set_constraints_binuclear import *
from get_binuclear_metalsite import *
from get_atoms_binuclear import *
from write_pdb_files import *

# Initialize object
no = get_binuclear_metalsite()
ga = get_atoms_binuclear()
scb = set_constraints_binuclear()
wpf = write_pdb_files()

# This is set for the constraint file writing.
SUBSTRATE_ATOMS_STR_0 = 'ZN1 O4 ZN2'
SUBSTRATE_ATOMS_STR_1 = 'ZN2 O4 ZN1'

# Cut-off distance used to define metal ligand residues
DISTANCE = 2.7


# Check it is working correct
def get_union_list(lst):
    unique_hetatom = ''
    for elem in range(len(lst) -1):
        a = lst.pop()
        if a in lst:
            unique_hetatom = a
    return unique_hetatom

def clean():
    # remove redundant files
    print 'clean'

def compute_cst(zn_coor,file_object,active_site,new_active_site):
    # Contains the block of constraint for each residue
    coordinates = []

    # Generate REMARK section for pdb file
    pdbfile =[]

    # Get the protein residues coordinating each metal atom
    # dummy variable to get the right metal ions
    j = 0

    # what is going on here?
    keylist = zn_coor.keys()
    keylist.sort()


    # For remarks section
    dummy = 1;      

    for key in keylist:
        ac_site,protein_lig = no.get_protein_ligand_single_metal(file_object,zn_coor[key],active_site,DISTANCE)

        # 15-08-2012 We use an array now
        decomp_zn_coor = zn_coor.values()
                
        if j == 0:

            # Values from the dictionary
            # Test that the first value corresponds to the metal ion

            if cross(decomp_zn_coor[0],zn_coor[key]).all() == array([0,0,0]).all():
                
                metal_1 = decomp_zn_coor[0]
                metal_2 = decomp_zn_coor[1]
            else:
                metal_2 = decomp_zn_coor[0]
                metal_1 = decomp_zn_coor[1]

            METAL = j
            
        elif j == 1:
        
            if cross(decomp_zn_coor[0],zn_coor[key]).all() == array([0,0,0]).all():
                metal_1,metal_2 = zn_coor.values()
            else:
                metal_2,metal_1 = zn_coor.values()
            METAL = j

        # What is going on here?

        tmp_constraints = []
        
        for i in protein_lig:
            info   = i.split()
            resn   = info[0]
            resid  = info[1]
            atom   = info[2]
            # Contains constraints prevents double residues
            t_c = str(resn+resid)
        
            protein_ligand_coordinates = scb.get_coordinates(file_object,resn,resid)
            
            dis,angA,angB,torA,torAB,torB = scb.get_geometry(resn,new_active_site,protein_ligand_coordinates,metal_1,metal_2)

            if j == 0:
                STRING = SUBSTRATE_ATOMS_STR_0
            else:
                STRING = SUBSTRATE_ATOMS_STR_1

            
            if float(dis) < DISTANCE and t_c not in tmp_constraints:
                coordinates.append(scb.write_constraint_file(resn,resid,atom,dis,angA,angB,torA,torAB,torB,STRING))
                pdbfile.append(scb.set_remarks_pdb(resn,resid,dummy))
                dummy = dummy +1

        j = j+1

    return coordinates, pdbfile



# Requires PDB file as input
# Returns constraint file, pdb-inputfile, 
def generate_output_rosetta(PDB,PDBNAME):

    # PDB file to be aligned to - here make loop over all files
    file_object = no.get_obj_pdbfile(PDB)

    # ZN coordinates and heteroatoms are retrived
    zn_coor,het = no.get_metal_ions(file_object)

    # Get keys of metal ions in dictionary
    # list with hetatoms
    tmp = []

    for item in zn_coor.items():
        # TMP contains all the hetero atoms
        tmp.extend(no.get_hetero_around_single_metal(item,het))

    # Here we determine which is the bridging water molecule
    # by getting the union between the two metal ions
    bridge_atom = get_union_list(tmp)
    


    # Returns heteroatom around Zn atom as dictionary
    active_site = no.get_hetero_around_metal(zn_coor,het,DISTANCE)

    if len(bridge_atom) == 0:
        
        # Due to a split later in the code
        # need to modify this as this could problems
        # This files is generated during the alignment - is this smart?
        read_ac = open('bridge_atom.pdb','r')

        for line in read_ac:
            ba = line
            
        # Getting vector with coordinates
        tmp_ba = ba.split()

        ba_coor = array([float(tmp_ba[6]),float(tmp_ba[7]),float(tmp_ba[8])])

        bridge_atom = 'nw\ bridge'

        active_site[bridge_atom] = ba_coor

    # Get the coordinates of the bridging oxygen
    new_active_site = active_site[bridge_atom]

    coordinates,pdbfile = compute_cst(zn_coor,file_object,active_site,new_active_site)

    # Write to file
    wpf.write_coordinates(coordinates,'constraints.cst')
    # Getting protein coordinates
    pdbfile.append(wpf.get_atoms_pdb(PDB))
    # Getting ligand
    pdbfile.append(wpf.get_ligand('chi_'+PDBNAME+'.pdb'))
    # Write pdbfile with ligand
    wpf.write_coordinates(pdbfile,'rosetta_'+PDBNAME+'.pdb')


def main():
    path = './'
    files = os.listdir(path)
    for fl in files:
        if os.path.isdir(fl):
            os.chdir(fl)
            PDB = fl+'.pdb'
            PDBNAME = str(fl)
            print 'doing PDB file',PDBNAME
            generate_output_rosetta(PDB,PDBNAME)
            print 'DONE',PDB
            os.chdir('../')


if __name__ == "__main__":
    main()

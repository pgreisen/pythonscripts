#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil
import argparse

'''
Fix_pdbfile get chain A as well as hetatms connected to chain A.
It takes all a pdb file as input and returns
the same file.

@parameter f: PDB file with metal site
@type       : string
@parameter m: Metal ion - default ZN 
@type       : string
@parameter n: Name of rosetta output pdb file 
@type       : string
@parameter l: Three ligand names with Zn and its first ligand
@type       : string
@parameter a: Coordinate file of aligned ligand and is append to the Rosetta input
@type       : string
@return     : constraint file, pdb file to use as input for rosetta


The default metal-protein ligand distance is set to 2.9 AA and for non-protein it is set to 2.7 AA


'''


class MetalSiteGeometry:

    def __init__(self):
        # Parameters for the distance geometry interactions
        # Distance between metal and protein
        self.DISTANCEMETAL = 2.9
        # Distance between metal and hetero atom
        self.DISTANCEHET = 2.7
        self.metal_coor = {}
        self.METAL = "ZN"
        self.coor_residues = ['THR','SER','LYS','ARG','ASN','HIS','GLN','TYR','ASP','GLU','CYS','MET','KCX']
        self.coor_atms = ['NE','NH2','NH1','NZ','OG','OG1','ND1','ND2', 'NE2','OE1', 'OE2','OD1','OD2','OH','SD','SG','O1','O2']
        self.protein_ligating_atoms = ['NE','NH2','NH1','NZ','OG','OG1','ND1','ND2', 'NE2','OE1', 'OE2','OD1','OD2','OH','SD','SG','O','O1','O2']
        self.coordination_sites_protein = {}
        self.chain = ""
        self.pdbfile = ""
        # LIGAND NAMES - three residues for the cst file
        self.ligandnames = ""
        self.ligandresname = ""
        self.metal_coordination = {}


    def write_stats_file(self):
        with open("stats.dat", "w") as f:
            f.write("Metal site,residue,distance,angle,torsion\n")
            for key in self.metal_coordination:
                for key2 in self.metal_coordination[key]:
                    if(len(self.metal_coordination[key][key2]) == 3):
                        f.write(str(key)+","+str(key2)+","+self.metal_coordination[key][key2][0]+","+self.metal_coordination[key][key2][1]+","+self.metal_coordination[key][key2][2]+"\n")
                    else:
                        f.write(str(key)+","+str(key2)+","+self.metal_coordination[key][key2][0]+"\n")



    # Requires PDB files
    # Returns list with lines in file
    # and sets the metal ions for later geometry determination
    def get_pdbfile(self):
        pdb = []
        with open(self.pdbfile,'r') as f:
            for line in f:
                if(line[0:4] == "ATOM"):
                    pdb.append( line )
                elif(line[0:4] == "HETA"):
                    atom_name = str(line[12:15]).strip()
                    het_id = line[13:14]
                    if atom_name == self.METAL:
                        # key is the residue name with the chain and residue number information
                        self.metal_coor[line[17:26]] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                        self.coordination_sites_protein[line[17:26]] = {}
                        self.metal_coordination[line[17:26]] = {}
                    pdb.append(line)
        return pdb

    # Requires list
    # Returns one list
    def merge(self,seq):
        merged = []
        for s in seq:
            for x in s:
                merged.append(x)
        return merged

    # Calculate angle between two vectors
    # from three points with vector2 being center
    # Return angle in degrees
    def get_calc_angle(self,vector1,vector2,vector3):
        x = vector1 - vector2
        y = vector3 - vector2
        dot_p = dot(x,y)
        norm_x = sqrt((x*x).sum())
        norm_y = sqrt((y*y).sum())
        # cosinus of angle between x and y
        cos_angle = dot_p / norm_x / norm_y
        angle = arccos(cos_angle)*57.3
        return angle

    # Return angle between two vectors in degrees
    def get_angle(self,n1,n2):
        dot_p = dot(n1,n2)
        norm_x = sqrt((n1*n1).sum())
        norm_y = sqrt((n2*n2).sum())
        cos_angle = dot_p / norm_x / norm_y # cosinus of angle between x and y
        angle = arccos(cos_angle)*57.3
        return angle

    # Requires 4 vectors (numpy.array([]))
    # Generates plane between four vectors
    # Return angle between the two planes
    def get_dihedral_angle(self,v1,v2,v3,v4):
        v12 = v2-v1
        v23 = v3-v2
        v34 = v4-v3
        # Normal vectors are generated
        n1 = cross(v12,v23)
        n2 = cross(v23,v34)

        angle = self.get_angle(n1,n2)
        # Getting the sign of the angle
        sign = dot(n1,v34)
        if sign < 0:
            angle = 360 - angle
        return angle

    # Requires pdb file
    # Returns protein atoms
    def get_atoms_pdb(self,PDB):
        atm = []
        for line in PDB:
            if line[0:4] == 'ATOM':
                atm.append(line)
        return atm

    # Requires pdb file
    # Returns heteroatoms
    def get_heteroatoms_pdb(self, PDB):
        atm = []
        for line in PDB:
            if line[0:4] == 'HETA':
                atm.append(line)
        return atm


    # Requires residue name of protein ligand
    # Returns names of residue required to determine
    # geometry (angle, torsion etc.)
    def get_list_of_atoms(self,resn):
        atms = {
            'HIS' : ['NE2','CE1','ND1'],
            'HIS2' : ['ND1','CE1','NE2'],
            'GLU' : ['OE1','CG','CD'],
            'GLU2' : ['OE2','CG','CD'],
            'ASP' : ['OD1', 'CG', 'CB'],
            'ASP2' : ['OD2', 'CG', 'CB'],
            'CYS' : ['SG','CB','CA'],
            'LYS' : ['NZ','CE','CD'],
            'KCX' : ['O1', 'C', 'NZ'],
            'THR' : ['OG1','CB','CA'],
            'SER' : ['OG','CB','CA'],
            'ARG' : ['NH1','CZ','NE'],
            'ARG2': ['NH2','CZ','NE'],
            'ARG3': ['NE','CZ','CD'],
            'GLN' : ['OE1','CD','CG'],
            'ASN' : ['OD1','CG','CB'],
            'TYR' : ['OH','CZ','CE2'],
            'MET' : ['SD','CG','CB']
        }
        return atms[resn]


    # Requires residue name of protein ligand
    # Returns names of residue required to determine
    # geometry (angle, torsion etc.)
    def get_list_of_atom_names(self,resn):
        atms = {
            'NE2' : 'NE2 CE1 ND1',
            'ND1' : 'ND1 CE1 NE2',
            'OE1' : 'OE1 CG CD',
            'OE2' : 'OE2 CG CD',
            'OD1' : 'OD1 CG CB',
            'OD2' : 'OD2 CG CB',
            'SG' :  'SG CB CA',
            'NZ' :  'NZ CE CD',
            'OG1' : 'OG1 CB CA',
            'OG' :  'OG CB CA',
            'NH1' : 'NH1 CZ NE',
            'NH2':  'NH2 CZ NE',
            'NE':  'NE CZ CD',
            'O1': 'O1 C NZ',
            'O2': 'O2 C NZ',
            #'GLN' : ['OE1','CD','CG'],
            #'ASN' : ['OD1','CG','CB'],
            'OH' : 'OH CZ CE2',
            'SD' : 'SD CG CB'
        }
        return atms[resn]


    # Get residue from pdb
    # returns dictionary with vector necessary for
    # generating constraints.
    def get_residue_constraint_pdb(self,PDB,resid, resn,chain):
        residue = {}
        # list of atoms of the protein
        try:
            atoms = self.get_list_of_atoms(resn)
            for line in PDB:
                res = line[17:20].rstrip()
                atm = line[0:4]
                if atm =='ATOM':
                    resnr = line[22:26].strip()
                    if res == resn[0:3]:
                        if resnr == resid:
                            atom = line[13:16].rstrip()
                            if atom in atoms and line[21:22].strip() == chain:
                                vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                                residue[atom] = vec
            return residue,atoms
        except:
            return "HOH", "O"

    # Requires list with pdb lines, metal coordinates ( set in the constructor )
    # Returns dictionary with ligands < DISTANCE from metal ion
    def get_protein_ligand_metal(self, PDB ):

        # loop over metal sites
        for metal_site in self.coordination_sites_protein.keys():
            TRUE = 0
            # Distance between ligands and metal ions
            mt_lig = []
            # coordinates of metal site
            metal_vec = self.metal_coor[metal_site]

            for line in PDB:

                res = line[17:20].rstrip()
                atm = line[0:4]

                if res in self.coor_residues and atm =='ATOM':
                    # Atom name of residue
                    atom = line[13:16].rstrip()
                    if atom in self.protein_ligating_atoms:
                        vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                        ds_lig = linalg.norm(vec-metal_vec)
                        if ds_lig < self.DISTANCEMETAL:
                            # assign value for ASP/GLU etc
                            if atom == 'OE2' and res == 'GLU':
                                res = 'GLU2'
                            elif atom == 'OD2' and res == 'ASP':
                                res = 'ASP2'
                            elif atom == 'ND1' and res == 'HIS':
                                res = 'HIS2'
                            iid = res+' '+line[21:22]+' '+line[22:26]
                            mt_lig.append(iid)
                            self.coordination_sites_protein[metal_site][iid] = vec


                elif atm == 'HETA':
                    vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                    ds_lig = linalg.norm(vec-metal_vec)
                    if ds_lig < self.DISTANCEMETAL and ds_lig > 0.001 :
                        iid = res+'  '+line[22:26]
                        mt_lig.append(iid)
                        self.coordination_sites_protein[metal_site][iid] = vec
                    else:
                        continue
                        print("Distance greater than cutoff:", ds_lig)

    # Map of amino acids necessary for
    # Rosetta internal naming
    # Requires residue name
    # Return rosetta 1-letter code
    def get_single_residue(self,resn):
        res_map = {
            'CYS' : 'C',
            'HIS' : 'H',
            'HIS2': 'H',
            'ASP' : 'D',
            'ASP2' : 'D',
            'GLU' : 'E',
            'GLU2' : 'E',
            'LYS' : 'K',
            'KCX' : 'K',
            'SER' : 'S',
            'ARG' : 'R',
            'ARG2' : 'R',
            'ARG3' : 'R',
            'THR' : 'T',
            'ASN' : 'N',
            'GLN' : 'Q',
            'TYR' : 'Y',
            'MET' : 'M'
            }
        return res_map[resn]


    # Map of atom names necessary for
    # Rosetta internal naming
    # Requires atom name
    # Return rosetta atom naming
    def get_single_atom(self,atom):
        Atm_map = {
            'SG' :  'S',
            'NE2' : 'Ntrp',
            'ND1' : 'Nhis',
            'OD1' : 'OOC',
            'OD2' : 'OOC',
            'OE1' : 'OOC',
            'OE2' : 'OOC',
            'NZ'  : 'Nlys',
            'OG1' : 'OH',
            'OG'  : 'OH',
            'NH1' : 'Narg',
            'NH2' : 'Narg',
            'NE'  : 'Narg',
            'OH'  : 'OH',
            'SD'  : 'S'
        }
        return Atm_map[atom]

    # Write constraint file with torsion A set to default value
    # Requires atom names, residue names, distance etc, ligand names
    # Returns string with parameters in Rosetta format
    def write_constraint_file(self,resn,resi,atom,chain,disAB, angB, torB, het_name_, chain_, residue_nr_,LIGANDRESNAME):
        atm_type     = self.get_list_of_atom_names( atom  )
        residue_type = self.get_single_residue(resn)

        const = '''

        # '''+het_name_+" "+chain_+" "+residue_nr_+''' - '''+str(resn)+' '+' '+chain+' '+str(resi)+'\n'+'''
        CST::BEGIN
        TEMPLATE::   ATOM_MAP: 1 atom_name:  '''+self.ligandnames+'''
        TEMPLATE::   ATOM_MAP: 1 residue3:  '''+LIGANDRESNAME+'''

        TEMPLATE::   ATOM_MAP: 2 atom_name: '''+str(atm_type)+''' ,
        TEMPLATE::   ATOM_MAP: 2 residue3:  '''+str(resn[0:3])+'''

        CONSTRAINT:: distanceAB:  '''+str(disAB)+'''  0.20  100.  1

        CONSTRAINT::    angle_B:  '''+str(angB)+'''   10.0  10.0  360.

        CONSTRAINT::  torsion_B:  '''+str(torB)+'''   10.0  10.00  360.

        CST::END'''
        
        return const

    # Computes geometry
    # Return gemetrical values
    def get_geometry(self,atoms, aminos, metal_site):

        if( len(atoms) == 3 ):
            l_a = atoms[0]
            l_s =atoms[1]
            l_t = atoms[2]
            dis = '%.2f' %(linalg.norm(self.metal_coor[metal_site] - aminos[l_a]))
            angB = '%.2f' %(self.get_calc_angle(self.metal_coor[metal_site], aminos[l_a], aminos[l_s]))
            torB = '%.2f' %(self.get_dihedral_angle(self.metal_coor[metal_site], aminos[l_a], aminos[l_s], aminos[l_t]))

            return dis, angB, torB
        else:
            dis =  '%.2f' %(linalg.norm( self.metal_coor[metal_site] - self.coordination_sites_protein[metal_site][aminos]  ))
            return dis

    # Requires PDB
    # Returns chain of first pdb
    def get_chain(self,pdb):
        for line in pdb:
            if line[0:4] == 'ATOM':
                chain = line[21:22]
                break
        return str(chain)


    # Requires residue name, type, number and chain
    # Returns string with constraint remarks
    def set_remarks_pdb(self,resn,resid,number,chain,metal_chain,metalname ):
        strng = 'REMARK   0 BONE TEMPLATE '+metal_chain+' '+metalname+'     0 MATCH MOTIF '+chain+' '+str(resn)+'    '+str(resid)+'   '+str(number)+'\n'
        return strng

    # Get geometrical constraint, write to file with them
    # Return remarks necessary for Rosetta PDB file
    def get_rosetta_constraint_files(self, PDB):
        remark = []
        dummy = 1
        rosetta_cst = []
        # get the primary interactions
        self.get_protein_ligand_metal( PDB )

        for metal_site in self.coordination_sites_protein.keys():
            for ligand in self.coordination_sites_protein[metal_site]:
                het_name_, chain_, residue_nr_ = metal_site.split()
                tm = ligand.split()
                resn =tm[0].rstrip()
                chain =tm[1].strip()
                resid =tm[2].rstrip()
                # there is no chain ID here
                aminos,atoms = self.get_residue_constraint_pdb(PDB,resid, resn,chain)

                # Define what the different values mean
                # Fix torsion
                if( len(atoms) == 3):
                    # returns what??
                    dis, angB, torB  = self.get_geometry(atoms, aminos, metal_site )
                    tmp_key = ligand+"_"+atoms[0]+"_"+atoms[1]+"_"+atoms[2]
                    self.metal_coordination[metal_site][tmp_key] = []
                    self.metal_coordination[metal_site][tmp_key].append(dis)
                    self.metal_coordination[metal_site][tmp_key].append(angB)
                    self.metal_coordination[metal_site][tmp_key].append(torB)


                    # reset metal site
                    LIGANDRESNAME = metal_site.split()[0]

                    # ToDo Add chain information in constraints and in remark
                    if( float(dis) <= self.DISTANCEMETAL):
                        # Adding modification
                        tmp_constraint = self.write_constraint_file( resn,resid,atoms[0],chain,dis,angB,torB,het_name_, chain_, residue_nr_,LIGANDRESNAME )
                        rosetta_cst.append(tmp_constraint) # write(tmp_constraint)

                        # Is this necessary still?
                        if resn == 'ARG2' or resn == 'ARG3':
                            resn = 'ARG'

                        elif resn == 'HIS2':
                            resn = 'HIS'

                        remark.append(self.set_remarks_pdb(resn,resid, dummy,chain,chain_, het_name_))
                        dummy = dummy +1
                else:
                    dis = self.get_geometry(atoms,ligand, metal_site)
                    self.metal_coordination[metal_site][ligand] = []
                    self.metal_coordination[metal_site][ligand].append(dis)

        return remark,rosetta_cst

    def main(self):

        parser = argparse.ArgumentParser(description="Generate restraints around metal site (Default: zinc ion")

        parser.add_argument('-m',dest='METAL',default='ZN',
                      help='Metal ion to get constraints for')
        parser.add_argument('-f',dest='pdbfile',
                      help='Name for output pdb file as input for Rosetta')
        parser.add_argument('-l',dest='ligandnames',default='ZN V1 V2',
                      help='Names from ligands involved in constraints e.g. ZN V1 V2')
        parser.add_argument('-t',dest='self.ligandresname',default='ZN',
                      help='Residue name of ligand default=ZN')
        parser.add_argument('-a',dest='LIGANDCOOR',default=False,
                      help='PDB coordinates of ligand appending it to rosetta pdb input')
        print("fine until here!")

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # PDB is the file object of the input file
        PDB = self.get_pdbfile()
        remark_pdb = open('rosetta_'+self.pdbfile,'w')
        remark, constraintfile = self.get_rosetta_constraint_files( PDB )

        tmp_pdb = self.get_atoms_pdb(PDB)
        pdb_file = remark + tmp_pdb

        with open("constraint.cst",'w') as f:
            for line in constraintfile:
                f.write(line)

        for line in pdb_file:
            remark_pdb.write(line)

        self.write_stats_file()


if __name__ == "__main__":
    run = MetalSiteGeometry()
    run.main()

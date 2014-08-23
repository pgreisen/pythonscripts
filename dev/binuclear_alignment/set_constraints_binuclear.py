#!/usr/bin/env python
import sys
from numpy import *

class set_constraints_binuclear:
    
    """
    Generate the data for the constraint file.
    1. Distances
    2. Angles
    3. Dihedral
    The torsion from the substrate is fixed into an idealized angle rest is taken from PDB
    
    Requires design file where 
    
    
    """

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
    # Is this used at all?
    def get_list_of_atoms(self,resn):        
        atms = {
            'SUB' : ['O2','C2', 'O4', 'O5', 'ZN1','O1','O6','06','O','O71','O3P','SG'],
            'HIS' : ['CB', 'CD', 'CE', 'CG','CE1', 'ND1', 'NE2'],
            'CYS' : ['CA', 'CB', 'SG'],
            'GLU' : ['CD', 'CG', 'OE1', 'OE2'],
            'ASP' : ['CB', 'CG', 'OD1', 'OD2'],
            'SER' : ['CA','CB','OG'],
            'GLN' : ['CD', 'CG', 'OE1'],
            'ASN' : ['CB', 'CG', 'OD1']
            
        }
        return atms[resn]


    def get_atm_type(self,atm_type):
        atm_tp = {
            'NE2' : 'NE2 CE1 ND1',
            'SG'  : 'SG CB CA',
            'OE1' : 'OE1 CD CG',
            'OE2' : 'OE2 CD CG',
            'OD1' : 'OD1 CG CB',
            'OD2' : 'OD2 CG CB',
            'ND1' : 'ND1 CG CB',
            'OG'  : 'OG CB CA'
            

            }
        return atm_tp[atm_type]
        



    def get_single_atom(self,atom):
        Atm_map = {
            'SG' : 'S',
            # Rosetta uses only one for the protonation
            #'NE2' : 'Ntrp',
            'NE2' : 'Nhis',
            'ND1' : 'Nhis',
            'OD1' : 'OOC',
            'OD2' : 'OOC',
            'OE1' : 'OOC',
            'OE2' : 'OOC',
            'OG'  : 'O'
        }

        return Atm_map[atom]

    def get_single_residue(self,resn):
        res_map = {
            'CYS' : 'C',
            'HIS' : 'H',
            'ASP' : 'D',
            'GLU' : 'E',
            'SER' : 'S',
            'GLN' : 'Q',
            'ASN' : 'N'
            }
        return res_map[resn]
         
    # Requires
    # Returns string with constraint remarks
    def set_remarks_pdb(self,resn,resid,number):
        strng = 'REMARK   0 BONE TEMPLATE X LG1    0 MATCH MOTIF A '+str(resn)+'    '+str(resid)+'   '+str(number)+'\n'
        return strng

    # Requires resn, resid, atom, distance, angleA, angleB. torAB, torB
    # Returns text constraint
    def write_constraint_file(self,resn,resi,atom,disAB,angA,angB,torA,torAB,torB,STRING):
        atm_type     = self.get_single_atom(atom)
        # Hack for ASN and GLN
        if resn == 'ASN':
            atm_type = 'ONH2'
        elif resn == 'GLN':
            atm_type = 'ONH2'

        atm_type_string = self.get_atm_type(atom)
        residue_type = self.get_single_residue(resn)
        const = '''

        # ZN - '''+str(resn)+' '+str(resi)+'\n'+'''
        CST::BEGIN
        TEMPLATE::   ATOM_MAP: 1 atom_name:  '''+STRING+'''
        TEMPLATE::   ATOM_MAP: 1 residue3:  LG1
        
        TEMPLATE::   ATOM_MAP: 2 atom_name: '''+str(atm_type_string)+''' 
        TEMPLATE::   ATOM_MAP: 2 residue3:  '''+str(resn)+'''
        
        CONSTRAINT:: distanceAB:  '''+str(disAB)+'''  0.20  200.  1
        CONSTRAINT::    angle_A:  '''+str(angA)+'''   25.00  10.0  180.
        CONSTRAINT::    angle_B:  '''+str(angB)+'''   10.0  10.0  180.
        
        CONSTRAINT::  torsion_A: '''+str(torA)+'''    10.0  1.00  360.
        CONSTRAINT:: torsion_AB: '''+str(torAB)+'''   10.0  2.00  360.
        CONSTRAINT::  torsion_B:  '''+str(torB)+'''   10.0  5.00  180.
        
        CST::END'''
        
        return const


    def get_atoms(self,resn,dicATM):
        return dicATM[resn]

    # Gets substrate and residue atomic coordinates
    # Requires file(which and what should it contain? right now I am using, resname, resid, substrate atoms(?), and atoms(?)
    # Returns coordinates
    # @Requires
    # @Returns
    #
    #    def get_coordinates(self,fname,resn,resid,sub_atms):
    # 16-08-2012
    # @Requires - filename, residue name and number
    # @Returns - dictionary with atom name and array of coordinates
    def get_coordinates(self,fname,resn,resid):

        atms = self.get_list_of_atoms(resn)
        rs = {};
        
        for line in fname:
            c_line = line[0:4]
            if c_line == 'ATOM':            
                tm = line.split()
                tmp_resid = str(line[23:26]).strip()
                tmp_resn = str(line[17:20])
                if tmp_resn == resn and tmp_resid == resid:
                    tmp_atom = str(line[13:16]).strip()
                    if tmp_atom in atms:
                        rs[tmp_atom] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])

        return rs



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
            'THR' : ['OG1','CB','CA'],
            'SER' : ['OG','CB','CA'],
            'ARG' : ['NH1','CZ','NE'],
            'ARG2': ['NH2','CZ','NE'],
            'ARG3': ['NE','CZ','CD'],
            'GLN' : ['OE1','CD','CG'],
            'ASN' : ['OD1','CG','CB'],
            'TYR' : ['OH','CZ','CE2']
            }
        return atms[resn]


    # Computes geometry
    # Return gemetrical values
    def get_geometry(self,resn,subcoor,aminos,metal_1,metal_2):

        if resn == 'HIS':
           if(linalg.norm(metal_1-aminos['ND1']) < linalg.norm(metal_1-aminos['NE2'])):
               atoms = self.get_list_of_atoms('HIS2')
           else:
               atoms = self.get_list_of_atoms('HIS')
        elif resn == 'ASP':
            if('OD1' in aminos.keys()):
                atoms = self.get_list_of_atoms('ASP')
            else:
                atoms = self.get_list_of_atoms('ASP2')               

        elif resn == 'GLU':
            if('OE1' in aminos.keys()):
                atoms = self.get_list_of_atoms('GLU')
            else:
                atoms = self.get_list_of_atoms('GLU2')               

        l_a = atoms[0]                                                                                  
        l_s =atoms[1]                                                                                   
        l_t = atoms[2]  

        dis =  '%.2f' %(linalg.norm(metal_1-aminos[l_a]))
        
        angA = '%.2f' %(self.get_calc_angle(subcoor,metal_1,aminos[l_a]))
        angB = '%.2f' %(self.get_calc_angle(metal_1,aminos[l_a],aminos[l_s]))
        
        torA = '%.2f'  %(self.get_dihedral_angle(metal_2,subcoor,metal_1,aminos[l_a]))

        torAB = '%.2f' %(self.get_dihedral_angle(subcoor,metal_1,aminos[l_a],aminos[l_s]))
        
        torB  = '%.2f' %(self.get_dihedral_angle(metal_1,aminos[l_a],aminos[l_s],aminos[l_t]))

        return dis,angA,angB,torA,torAB,torB

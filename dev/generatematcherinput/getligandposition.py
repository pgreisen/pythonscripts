#!/usr/bin/env python
import sys,os
from numpy import *


class GetLigandPositions:

    """
    Generate resfile for design for all residues XX Aangstrom
    away from ligand - default is 8
       
    @requires:
    
    python get_positions_around_ligand.py PDBFILE POSITIONFILE LIGANDATOMNAME
    
    @return:
    the positions with the cutoff
    
    """
    
    
    # DISTANCE_AWAY_FROM_LIGAND = 6


    def __init__(self):
        self.DISTANCE_AWAY_FROM_LIGAND = 6
    
    # @require a position file
    # returns a list with potential positions
    def get_positions(self,posfile):
        positions = []
        tmp = open(posfile,'r')
        for line in tmp:
            tmp2 = line.split()
            for i in tmp2:
                positions.append(i)
        tmp.close()
        return positions


    # @require a position file
    # returns a list with potential positions
    def get_pdbfile(self,pdbfile):
        lines = []
        tmp = open(pdbfile,'r')
        for line in tmp:
            lines.append(line)
        tmp.close()
        return lines


    # @require atom names of ligand
    # return an average position of the ligand atom
    # good with aromatic rings
    def get_ligandatom_aromatic_coordinates(self,pdbfile,names):
        ligand_coor = array([0,0,0])
        for line in pdbfile:
            if line[0:4] == 'HETA' and line[12:16].strip() in names:
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                ligand_coor = ligand_coor + array([float(x),float(y),float(z)])
        return ligand_coor/2


    # @require atom name of ligand
    # return ligand coordinates
    def get_ligandatom_coordinates(self,pdbfile,name):
        ligand_coor = []
        for line in pdbfile:
            if line[0:4] == 'HETA' and line[12:16].strip() == name:
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                ligand_coor.append(array([float(x),float(y),float(z)]))
        return ligand_coor

    # @require pdb-file, array with positions on protein scaffold, array with ligand coordinates
    # return coordinates of the positions
    def get_protein_positions(self,pdbfile,ligand_coor,positions):
        new_positions = []
        tmp_positions = []
        for line in pdbfile:
            tmp_pos = str(line[23:26]).strip()
            if line[0:4] == 'ATOM' and  tmp_pos in positions:
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                tmp_positions.append(array([float(x),float(y),float(z)]))

                if(linalg.norm(ligand_coor - array([float(x),float(y),float(z)])) <= self.DISTANCE_AWAY_FROM_LIGAND) and tmp_pos not in new_positions:
                    new_positions.append(tmp_pos)
    
        return new_positions

    # @require pdb-file, array with positions on protein scaffold, array with ligand coordinates
    # return coordinates of the positions
    def get_protein_positions(self,pdbfile,ligand_coor):
        new_positions = []
        tmp_positions = []
        for line in pdbfile:
            tmp_pos = str(line[23:26]).strip()
            if line[0:4] == 'ATOM':
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                tmp_positions.append(array([float(x),float(y),float(z)]))

                if(linalg.norm(ligand_coor - array([float(x),float(y),float(z)])) <= self.DISTANCE_AWAY_FROM_LIGAND):
                    new_positions.append(int(tmp_pos))

        return sorted(set(new_positions))


    # Get direction between ligand atom and protein residue - are they point towards
    # each other.
    def get_direction_residue_ligand_atom(calpha_coordinates, cbeta_coordinates, ligand_atom):
        dir = False
        vector1 = cbeta_coordinates - calpha_coordinates 
        vector2 = cbeta_coordinates - ligand_atom
        direction = dot(vector1,vector2)
        if direction < 0:
            dir = True
        return dir

#!/usr/bin/env python
import numpy
import math
from numpy.ma import arccos, dot


class GetLigandPositions:
    """
    Generate resfile for design for all residues XX Aangstrom
    away from ligand - default is 8
       
    @requires:
    
    python get_positions_around_ligand.py PDBFILE POSITIONFILE LIGANDATOMNAME
    
    @return:
    the positions with the cutoff
    
    """

    def __init__(self):
        self.DISTANCE_AWAY_FROM_LIGAND = 6

    # @require a position file
    # returns a list with potential positions
    #def get_positions(self, posfile):
    #    """
    #
    #    @param posfile:
    #    @return:
    #    """
    #    positions = []
    #    tmp = open(posfile, 'r')
    #    for line in tmp:
    #        tmp2 = line.split()
    #        for i in tmp2:
    #            positions.append(i)
    #    tmp.close()
    #    return positions

    # @requires two vectors
    # @returns angle between two vectors
    def get_angle_between_vectors(self, vector1, vector2):
        """
        Calculates the angle between two vectors
        @param vector1:
        @param vector2:
        """
        dot_p = dot(vector1, vector2)
        norm_x = math.sqrt((vector1 * vector1).sum())
        norm_y = math.sqrt((vector2 * vector2).sum())
        cos_angle = dot_p / norm_x / norm_y # cosinus of angle between x and y
        angle = arccos(cos_angle) * 57.3
        return angle


    def get_vector_atom_name(self, pdbfile, residue_number, residue_name, atom_name="CA"):

        """

        @param pdbfile: list with protein ligand coordinates
        @param residue_number:number of residue
        @param residue_name: name of residue
        @param atom_name:String of atom name

        """

        for line in pdbfile:
            if ( str(line[23:26]).strip() == residue_number):
                if ( str(line[17:20]).strip() == residue_name ):
                    if ( str(line[13:16]).strip() == atom_name):
                        x = numpy.str(line[30:38]).rstrip()
                        y = numpy.str(line[38:46]).rstrip()
                        z = numpy.str(line[46:54]).rstrip()
                        return numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)])

        pass


    def get_angle_between_ligand_and_protein(self, vector1, vector2, vector3):

        """

        @param vector1: c-alpha coordinates
        @param vector2: c-beta coordinates
        @param vector3: ligand atom
        @return: angle between ca-cb cb-la
        """
        x = vector1 - vector2
        y = vector3 - vector2
        # fix this at some point
        y = y[0]

        dot_p = numpy.dot(x, y)
        norm_x = math.sqrt((x * x).sum())
        norm_y = math.sqrt((y * y).sum())
        # cosinus of angle between x and y
        cos_angle = dot_p / norm_x / norm_y
        angle = arccos(cos_angle) * 57.3
        return angle


    def get_pdbfile(self, pdbfile):
        lines = []
        tmp = open(pdbfile, 'r')
        for line in tmp:
            lines.append(line)
        tmp.close()
        return lines


    # @require atom names of ligand
    # return an average position of the ligand atom
    # good with aromatic rings
    def get_ligandatom_aromatic_coordinates(self, pdbfile, names):
        ligand_coor = numpy.array([0, 0, 0])
        debug = []
        for line in pdbfile:
            if line[0:4] == 'HETA' and line[12:16].strip() in names:
                x = numpy.str(line[30:38]).rstrip()
                y = numpy.str(line[38:46]).rstrip()
                z = numpy.str(line[46:54]).rstrip()
                debug.append( numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)]) )
                ligand_coor = ligand_coor + numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)])
        print "aromatic center is determined to ", ligand_coor / len(names)
        return  ligand_coor / len(names)


    # @require atom name of ligand
    # return ligand coordinates
    def get_ligandatom_coordinates(self, pdbfile): #, name):
        ligand_coor = []
        for line in pdbfile:
            if line[0:4] == 'HETA':# and line[12:16].strip() == name:
                x = numpy.str(line[30:38]).rstrip()
                y = numpy.str(line[38:46]).rstrip()
                z = numpy.str(line[46:54]).rstrip()
                ligand_coor.append(numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)]))
        return ligand_coor

    # @require pdb-file, array with positions on protein scaffold, array with ligand coordinates
    # return coordinates of the positions
    def get_positions(self, pdbfile, ligand_coor):
        """
        @rtype : list
        @param pdbfile: list with pdbfile
        @param ligand_coor: the atom coordinates of the ligand atom
        @param positions: the positions identify
        @return: the positions where the Ca-Cb Cb-LigandAtom have an angle >= 60 and <= -60
        """
        new_positions = []
        tmp_positions = []

        for line in pdbfile:
            tmp_pos = str(line[23:26]).strip()
            if( line[0:4] == 'ATOM'):
                x = numpy.str(line[30:38]).rstrip()
                y = numpy.str(line[38:46]).rstrip()
                z = numpy.str(line[46:54]).rstrip()
                # get type of residue
                residue_name = str(line[17:20])
                residue_number = str(line[23:26]).strip()
                #tmp_positions.append(numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)]))
                # loop over the coordinates from the ligand
                for j in ligand_coor:
                    if (numpy.linalg.norm(j - numpy.array([numpy.float(x), numpy.float(y), numpy.float(z)])) <= self.DISTANCE_AWAY_FROM_LIGAND) and tmp_pos not in new_positions:
                        # If residue is not equal to glycine we need to determine the angle
                        if (residue_name != "GLY"):
                            # Get the c-alpha and the c-beta vector
                            calpha_vector = self.get_vector_atom_name(pdbfile, residue_number, residue_name, "CA")
                            cbeta_vector = self.get_vector_atom_name(pdbfile, residue_number, residue_name, "CB")
                            direction_angle = self.get_angle_between_ligand_and_protein(calpha_vector, cbeta_vector, ligand_coor)
                            # print "The directional angle is ", direction_angle
                            if ( 60.0 <= direction_angle ):
                                new_positions.append(tmp_pos)
                            elif (  direction_angle <= -60.0 ):
                                new_positions.append(tmp_pos)
                        elif (residue_name == "GLY" ):
                            new_positions.append(tmp_pos)

        return new_positions
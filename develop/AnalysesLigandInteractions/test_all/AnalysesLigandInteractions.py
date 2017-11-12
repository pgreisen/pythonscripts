#!/usr/bin/env python
import sys
from numpy import *
import argparse
from scipy import linalg

"""
Generate resfile for enzyme design

@requires pdb file with ligand present
@return resfile with natro for all residue > 10 AA and nataa > 8 AA

"""

class GenerateResfileAroundResidues:


    def __init__(self):
        self.pdbname = ""

        self.DISTANCE1 = 6
        self.DISTANCE2 = 8
        self.residuenumber = ""
        self.residues = []
        self.pdbfile = []
        self.datafile = []
        self.carbohydrates = ["BGC", "GLC"]
        self.carbohydrates_coordinates = {}
        self.carbohydrate_links = []
        self.max_dis = 4.5
        self.max_polar_interaction = 3.5
        self.polar = ['O','N']

        self.polar_geometry = []
        self.packing = []
        self.hydrogen_bond_cutoff = 2

        self.carb_output = {}


    def get_pdbfile(self,file):
        with open(file,'r') as f:
            for line in f:
                if( line[0:4] == "ATOM" ):
                    self.pdbfile.append(line)

                if( line[0:4] == "HETA" ):
                    # test residue is carbohydrate of interest
                    if line[17:20].strip() in self.carbohydrates:
                        test_key = line[17:20]+"_"+line[21:22]+"_"+line[22:26].strip()
                        if test_key not in self.carbohydrates_coordinates.keys():
                            self.carbohydrates_coordinates[test_key] = {}
                        atomname = line[13:15]
                        x,y,z = self.get_xyz(line)
                        self.carbohydrates_coordinates[test_key][atomname] = array([x,y,z])
                        self.pdbfile.append(line)

                if( line[0:4]  == "LINK" ):


                    for carb in self.carbohydrates:

                        if( line.find(carb) != -1 ):

                            self.carbohydrate_links.append( line )

        # print self.carbohydrates_coordinates


    def get_xyz(self,line):
        x = str(line[30:38]).rstrip()
        y = str(line[38:46]).rstrip()
        z = str(line[46:54]).rstrip()
        return float(x), float(y), float(z)


    # return coordinates of residue
    # get_residue_coordinates
    def get_residue_coordinates(self):
        # Debug
        ligand_coor = []
        for line in self.pdbfile:
            if line[22:26].strip() in self.residues:
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                ligand_coor.append(array([float(x),float(y),float(z)]))
        return ligand_coor

    def get_length_protein(self):
        first = '0'
        start = ''
        end = ''
        for line in self.pdbfile:
            if line[0:4] == 'ATOM':
                if first == '0':
                    start = str(line[23:26]).rstrip()
                    first = '1'
                elif first == '1':
                    end = str(line[23:26]).rstrip()
        return int(start),int(end)


    def generate_geometry_analysis(self):
        # loop over sugar groups
        for carb_unit in self.carbohydrates_coordinates.keys():
            for atom in self.carbohydrates_coordinates[carb_unit]:
                # loop over pdb coordinates in pdbfile
                for line in self.pdbfile:
                    if( line[0:4] == "ATOM"):

                        x, y, z = self.get_xyz(line)
                        prt_vector = array([x, y, z])
                        dis = linalg.norm( prt_vector - self.carbohydrates_coordinates[carb_unit][atom] )

                        # make key and take length of list

                        if( dis <= self.max_dis ):

                            tmp_atom_name = line[13:17].strip()
                            tmp_resname = line[17:20].strip()
                            tmp_resnr = line[22:26].strip()

                            new_key = self.pdbname+'_'+carb_unit
                            # print new_key
                            if( new_key not in self.carb_output.keys() ):
                                self.carb_output[new_key] = {}
                                self.carb_output[new_key]["polar"] = []
                                self.carb_output[new_key]["pack"] = []



                            # print carb_unit, atom, dis, tmp_atom_name, tmp_resname, tmp_resnr
                            # get polar interactions
                            if( atom[0] in self.polar and tmp_atom_name[0] in self.polar and dis <= self.max_polar_interaction  ):
                                tmp_polar_str = carb_unit+','+atom+','+str(round(dis,2))+','+tmp_atom_name+','+tmp_resname+','+tmp_resnr

                                self.carb_output[new_key]["polar"].append(tmp_polar_str)

                                self.polar_geometry.append( tmp_polar_str )
                                # print carb_unit, atom, dis, tmp_atom_name, tmp_resname, tmp_resnr
                            else:
                                tmp_pack_str = carb_unit + ',' + atom + ',' + str( round(dis, 2)) + ',' + tmp_atom_name + ',' + tmp_resname + ',' + tmp_resnr
                                self.packing.append( tmp_pack_str )
                                self.carb_output[new_key]["pack"].append(tmp_pack_str)





    def main(self):
        parser = argparse.ArgumentParser(description="Generate Residue File Around A List of Residues")
        parser.add_argument('-f',dest='datafile', help='PDB files' )
        parser.add_argument('--residuenumber',dest='residuenumber', help='Residues with PDB numbering which are used to generate the resfile' )
        parser.add_argument('--distance1',dest='distance1', default=6,type=int, help="If atoms are below this value the residue is allowed to repack")
        parser.add_argument('--distance2',dest='distance2', default=8,type=int, help="Above this threshold residues are not allowed to repack")
        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # Get pdbfile
        #self.get_pdbfile( self.datafile)

        self.get_pdbfile( self.datafile )

        self.pdbname = self.datafile.split('.')[0]
        ##print self.carbohydrate_links

        ##print self.carbohydrates_coordinates

        self.generate_geometry_analysis()

        pol_out = []
        pack_out = []

        for keykey in self.carb_output.keys():
            # print keykey, self.carb_output[keykey]["polar"]

            if( len( self.carb_output[keykey]["polar"] ) >= self.hydrogen_bond_cutoff ):
                for i in self.carb_output[keykey]["polar"]:
                    pol_out.append(i)
                for j in self.carb_output[keykey]["pack"]:
                    pack_out.append(j)


        with open(self.pdbname+"_polar.csv",'w' ) as f:
            for line in pol_out:
                f.write(line+'\n')

        with open(self.pdbname + "_pack.csv", 'w') as f:
            for line in pack_out:
                f.write(line + '\n')


        
if __name__ == "__main__":
    run = GenerateResfileAroundResidues()
    run.main()

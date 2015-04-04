import os,sys, argparse
from numpy import *
from numpy import linalg as LA
from collections import defaultdict


class GeneratePositionFile:


    def __init__(self):
        self.one_chain_only = False
        self.angle  = 0.0
        self.dis = 2.0
        self.maximum_distance = 9
        self.backbone_coordinates = {}
        self.cbeta_dummy = {}
        self.positions = []
        self.general_format = 0
        self.chainA = {}
        self.chainB = {}
        self.offset = 0
        self.posA = []
        self.posB = []




    def get_cbeta_dummy_vector(self, n_atom, ca_atom, c_atom):
        # interface_vector_calculate.cc l. 460
        '''
        // 		//figure out where CB would be
        // 		xyzVector< Real > v1( midpoint( C_xyz, N_xyz)  );
        // 		xyzVector< Real > v2( 2*CA_xzy - v1 ); //reflection of v1 about CA position
        // 		Real d1( magnitdue( v1 - C_xyz) );  //distance from midpoint (v1 or v2) to an atom
        // 		xyzVector< Real > dir_CB( cross(N_xyz - CA_xyz, C_xyz - CA_xyz).normalize() ); //direction of CB from V2
        // 		xyzVector< Real > CB_xyz( v2 + dir_cb * d1);
        // 		return CB_xyz;
        '''

        vec1 = (c_atom + n_atom)/ 2
        vec2 = (2*ca_atom - vec1)
        # 2norm of vector
        d1 = LA.norm( vec1 - c_atom )
        dir_cb = cross(n_atom - ca_atom, c_atom -ca_atom )
        d_cb = dir_cb / LA.norm( dir_cb )
        cb_displace  = vec2 + dir_cb * self.dis;
        return cb_displace

    def get_backbone_coordinates(self,pdbfile):
        backbone = [ "N", "CA", "C" ]
        with open(pdbfile, 'r') as f:
            for line in f:
                if( len(line) > 4 and line[0:4] == "ATOM" ):

                    tmp = line.split()
                    key = tmp[3]+"_"+tmp[4]+"_"+tmp[5]

                    # debug
                    # print key # output (residue name, chain, residue number ): ALA_A_200

                    # initialize key if it not already present in the dictionary
                    if( self.backbone_coordinates.has_key( key ) == False):
                        self.backbone_coordinates[key] = {}

                    # If not a glycine use the already assigned c-beta atom coordinates
                    if( line[17:20] != "GLY" and tmp[2] == "CB" ):
                        self.cbeta_dummy[key] = array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])
                        # debug
                        # print array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])


                    elif( line[17:20] == "GLY" ):
                        self.backbone_coordinates[key][tmp[2]] = array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])

                    # we do not care about the rest of the atoms in the file
                    else:
                        continue
                        #print line
                        #sys.exit("Strange behavior for generating position file")



    def get_backbone_coordinates_between_chain_A_B(self,pdbfile):

        backbone = [ "N", "CA", "C" ]
        with open(pdbfile, 'r') as f:
            for line in f:

                if( len(line) > 4 and line[13:15] == "CA" and str(line[21:22]) == 'A' ):
                    tmp = line.split()
                    key = tmp[3]+"_"+tmp[4]+"_"+tmp[5]
                    # initialize key if it not already present in the dictionary
                    if( self.chainA.has_key( key ) == False):
                        self.chainA[key] = {}

                    self.chainA[key] = array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])

                elif( len(line) > 4 and line[13:15] == "CA" and str(line[21:22]) == 'B' ):
                    tmp = line.split()
                    key = tmp[3]+"_"+tmp[4]+"_"+tmp[5]

                    # initialize key if it not already present in the dictionary
                    if( self.chainB.has_key( key ) == False):
                        self.chainB[key] = {}

                    self.chainB[key] = array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])



    def get_positions(self):
        '''
        self.posA = []
        self.posB = []
        '''
        for key in self.chainA:
            for key_two in self.chainB:
                # print  self.chainA[key], self.chainB[key_two]
                dis = LA.norm( self.chainA[key] - self.chainB[key_two])
                if( dis <= self.maximum_distance):

                    self.posA.append( key.split('_')[2])
                    self.posB.append( key_two.split('_')[2])


    def write_position_file(self,posA, posB,name="AB.pos"):
        with open( name ,'w') as f:
            f.write("N_CST 2\n 1 : ")
            for i in set( posA ):
                f.write(str( int(i) + self.offset)+' ')
            f.write("\n 2 : ")
            for j in set( posB ):
                f.write(str( int(j) + self.offset)+' ')


    def write_general_position_file(self,posA,posB,name="pos.pos"):

        with open( name ,'w') as f:
            for i in set( posA ):
                f.write(str( int(i) + self.offset)+' ')
            for j in set( posB ):
                f.write(str( int(j) + self.offset)+' ')


    def main(self):

        parser = argparse.ArgumentParser(description="Generate a position file to be used in connection with matching")
        # get the initial rosetta design as input
        parser.add_argument("-d","--distance", dest="distance", help="Maximum distance to be included in the search for Cbeta-Cbeta distance (Default 5 AA)", type=float)

        parser.add_argument("-b", "--bundles", dest="helical_bundle", help="Four chains helical bundle with four chains is set to true", action="store_true", default=False )

        parser.add_argument("-f", "--pdbfile", dest="pdbfile", help="The name of the pdb file", default=None, type=str )

        parser.add_argument("--offset", dest="offset", help="Offset to be added to the position", default=0, type=int )

        parser.add_argument("--format", dest="general_format", help="Matcher format is default but none zero value you will get the general format", default=0, type=int )

        input_variables = parser.parse_args()

        # old version
        # changed 04-12-2014
        # pdbfile = sys.argv[1]

        pdbfile = input_variables.pdbfile
        if( input_variables.distance):
            self.maximum_distance = input_variables.distance

        if( input_variables.helical_bundle):
            self.one_chain_only = input_variables.helical_bundle


        self.offset = input_variables.offset

        self.get_backbone_coordinates_between_chain_A_B( pdbfile )
        self.get_positions()
        if(input_variables.general_format == 0):

            self.write_position_file(self.posA, self.posB, "AB.pos")
            self.write_position_file(self.posB, self.posA, "BA.pos")

        else:
            self.write_general_position_file(self.posA, self.posB, "pos.pos")


if __name__ == "__main__":
    run = GeneratePositionFile()
    run.main()



"""
   def test_script(self):

        ATOM    199  N   ALA A  13       2.193  -4.925  94.792  1.00  0.00
        ATOM    200  CA  ALA A  13       1.423  -4.056  93.919  1.00  0.00
        ATOM    201  C   ALA A  13       2.277  -3.489  92.794  1.00  0.00
        ATOM    202  O   ALA A  13       1.874  -3.470  91.642  1.00  0.00
        ATOM    203  CB  ALA A  13       0.791  -2.932  94.724  1.00  0.00
        '''

        n_atom_13 = array([2.193,  -4.925,  94.792])
        ca_atom_13 = array([1.423,  -4.056,  93.919])
        c_atom_13 = array([2.277,  -3.489,  92.794])
        cb_atom_13 = array([0.791, -2.932, 94.724])

        nv = self.get_cbeta_dummy_vector(n_atom_13, ca_atom_13, c_atom_13 )

        ATOM    233  N   GLU A  15       5.087  -4.737  91.604  1.00  0.00
        ATOM    234  CA  GLU A  15       5.405  -5.851  90.722  1.00  0.00
        ATOM    235  C   GLU A  15       4.279  -6.124  89.730  1.00  0.00

        n_atom_15 = array([5.087, -4.737, 91.604 ])
        ca_atom_15 = array([5.405, -5.85, 90.722])
        c_atom_15 = array([4.279, -6.124, 89.730])

        n15 = self.get_cbeta_dummy_vector(n_atom_15, ca_atom_15, c_atom_15 )

        print "n15:  ",n15

        # the values are taken from 3helix_1.pdb
        # residue 17, 18, 19
        ##n_atom = array([1.980, -3.951, 88.768])
        ##ca_atom = array([1.969, -2.814, 87.845])
        ##c_atom = array([3.010, -2.970, 86.739])
        ##cb_atom = array ([2.224, -1.481, 88.590])

        ##nd = cb_atom + LA.norm(cb_atom - ca_atom) #* self.dis

        ##nv = self.get_cbeta_dummy_vector(n_atom, ca_atom, c_atom )
"""

import os,sys
from numpy import *
from numpy import linalg as LA


class GeneratePositionFile:


    def __init__(self):
        self.one_chain_only = False
        self.angle  = 0.0
        self.dis = 2.0
        self.maximum_distance = 5
        self.backbone_coordinates = {}
        self.cbeta_dummy = {}
        self.positions = []


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



    def get_cbeta_position(self):
        '''
        Method to generate a pseudo atom for the c-beta atom. Only used for glycine residues
        '''

        for key in self.backbone_coordinates:
            if( str( key[0:3] ) == "GLY"):
                self.cbeta_dummy[key] = self.get_cbeta_dummy_vector( self.backbone_coordinates[key]["N"], self.backbone_coordinates[key]["CA"], self.backbone_coordinates[key]["C"]  )
                print self.cbeta_dummy[key]



    def get_positions(self):

        for key in self.cbeta_dummy:
            for key_two in self.cbeta_dummy:

                # remove identical keys
                if( key == key_two):

                    continue
                # remove residue in the same chain
                elif (key.split('_')[1] == key_two.split('_')[1] and self.one_chain_only):
                    continue
                else:
                    dis = LA.norm(self.cbeta_dummy[key] - self.cbeta_dummy[key_two])
                    if( dis <= self.maximum_distance):
                        self.positions.append( key.split('_')[2])
                        break


    def main(self):
        pdbfile = sys.argv[1]
        self.get_backbone_coordinates(pdbfile)
        self.get_cbeta_position()
        self.get_positions()

        print self.backbone_coordinates
        print "cbeta positions",self.cbeta_dummy
        print "positions",self.positions


        with open("new_position.pos", 'w') as f:
            for line in self.positions:
                f.write(line+' ')


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

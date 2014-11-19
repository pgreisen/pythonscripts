import os,sys

from numpy import *
from numpy import linalg as LA


class GeneratePositionFile:


    def __init__(self):
        self.angle  = 0.0
        self.dis = 2.0
        self.maximum_distance = 7
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
                tmp = line.split()
                key = tmp[3]+"_"+tmp[4]+"_"+tmp[5]
                if( self.backbone_coordinates.has_key( key ) == False):
                    self.backbone_coordinates[key] = {}
                if( tmp[2] in backbone ):
                    self.backbone_coordinates[key][tmp[2]] = array([ float(tmp[6]), float(tmp[7]), float(tmp[8])])


    def get_cbeta_position(self):

        for key in self.backbone_coordinates:
            self.cbeta_dummy[key] = self.get_cbeta_dummy_vector( self.backbone_coordinates[key]["N"], self.backbone_coordinates[key]["CA"], self.backbone_coordinates[key]["C"]  )


    def get_positions(self):

        for key in self.cbeta_dummy:

            for key_two in self.cbeta_dummy:
                # remove identical keys
                if( key == key_two):
                    continue
                # remove residue in the same chain

                elif (key.split('_')[1] == key_two.split('_')[1]):
                    # print key.split('_')[1], key_two.split('_')[1]
                    continue
                else:
                    dis = LA.norm(self.cbeta_dummy[key] - self.cbeta_dummy[key_two])
                    if( dis <= self.maximum_distance):
                        self.positions.append( key.split('_')[2])
                        break





    def main(self):

        pdbfile = sys.argv[1]

        self.get_backbone_coordinates(pdbfile)

        #
        self.get_cbeta_position()
        self. get_positions()

        with open("position.pos", 'w') as f:
            for line in self.positions:
                f.write(line+'+')


        # print "The number of residues in the position file is: ", len(self.positions)


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
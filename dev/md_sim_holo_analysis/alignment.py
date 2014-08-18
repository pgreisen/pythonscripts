from numpy import *
import numpy

'''
Superposition of two matrices using the algorithm
develpoed by Kabsch and Wolfgang


Make an alignment of C-alpha atoms only

'''

class Alignment:

    def get_data(self,pdbfile):
        # Coordinates
        dt = []
        # atom names
        a_n = []
        # Looping over line in file
        for line in pdbfile:
            sl = line.split()
            x = float(line[31:38])
            y = float(line[38:46])
            z = float(line[46:54])
            tmp = x,y,z  
            dt.append(tmp)
            a_n.append(sl[2])

        return array(dt),a_n


    def get_calpha_coordinates(self,pdbfile):
        # Coordinates
        """

        @param pdbfile:
        @return:
        """
        dt = []
        # Looping over line in file
        for line in pdbfile:
            sl = line.split()
            # Only take c-alpha coordinates
            if( line[13:15] == "CA"):
                x = float(line[31:38])
                y = float(line[38:46])
                z = float(line[46:54])
                tmp = x,y,z
                dt.append(tmp)

        return array(dt)


    def write_pdb(self,atomnames,coordinates,fname='superimposed.pdb'):
        '''Write a pdb file of the aligned coordinates '''
        al = len(atomnames)
        # Name of chain
        ch = '  SUB'
        wname = open(fname,'w')
        for i in range(al):
            wname.write('ATOM   '+str(i).rjust(3)+'   '+atomnames[i]+' '+ch+' X '+str(i).rjust(3)+'    '+str(coordinates[i][0])[0:6]+'  '+str(coordinates[i][1])[0:6]+'  '+str(coordinates[i][2])[0:6]+'\n')
            

    # Kabsch Algorithm for superimposing
    # Superimpose matrix1 on matrix2
    def get_rotate_translate(self,matrix1,matrix2):

        assert shape(matrix1) == shape(matrix2), 'Matrices not of same dimensions'
        
        # Store number of rows
        nrows = shape(matrix1)[0]
        
        # Getting centroid position for each selection
        avg_pos1 = matrix1.sum(axis=0)/nrows
        avg_pos2 = matrix2.sum(axis=0)/nrows


        # Translation of matrices
        avg_matrix1 = matrix1-avg_pos1
        avg_matrix2 = matrix2-avg_pos2

        # Covariance matrix
        covar = numpy.dot(numpy.transpose(avg_matrix1), avg_matrix2)


        # Do the SVD in order to get rotation matrix
        u,s,wt = linalg.svd(covar)
        
        # Rotation matrix
        # Transposition of u,wt
        rot_matrix = wt.T*u.T
        
        # Insure a right-handed coordinate system
        if linalg.det(rot_matrix) > 0:
            wt[2] = -wt[2]
            rot_matrix = transpose(dot(transpose(wt),transpose(u)))

        trans_matrix = avg_pos2-dot(avg_pos1,rot_matrix)
        return trans_matrix, rot_matrix


    # Returning the superimposed coordinates
    # Requires two arrays same length
    # Returns array of transformed coordinates
    # A is superimposed on B
    def get_transformed_coor(self,a,b):
        tr, rt = self.get_rotate_translate(a,b)
        nw_coor = dot(a,rt) + tr
        return nw_coor

    def transform_coordination(self,proteinA, t_m, r_m):
        nw_coor = dot(proteinA , r_m) + t_m
        return nw_coor


    # Executes the above scripts
    def align_proteinA_onto_proteinB(self,proteinA, proteinB):

        pa =  self.get_calpha_coordinates( proteinA )
        pb =  self.get_calpha_coordinates( proteinB )

        t_m , r_m = self.get_rotate_translate( pa, pb )

        pa_new_coor = self.get_transformed_coor(pa, pb)

        print pa_new_coor



    def rmsd(self, V, W):
        # coordinates
        D = 3 # len(V[0][0])
        N = len( V )
        rmsd = 0.0

        # for v, w in zip(V, W):
        for key in V:

            rmsd += sum([(V[key]-W[key])**2.0 for i in range(D)])
            #rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])

        return numpy.sqrt(rmsd/N)



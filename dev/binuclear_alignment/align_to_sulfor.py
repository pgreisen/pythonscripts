from numpy import *

''' Superimpose one coordinate set on another
    The algorithm is based on Landau mechanics
    Returns the coordinates in pdb file
'''

class align_to_substrate:
    # Transform coordinates into numpy array
    # Requires coordinates in pdb format
    # Returns numpy array with coordinates
    def get_data(self,filename):
        # Need to add try/catch
        fname = open(filename,'r')
        # Coordinates
        dt = []
        # atom names
        a_n = []
        # Looping over line in file
        for line in fname:
            sl = line.split()
            # Bug fix 28-12-2009
            x = float(line[31:39])
            y = float(line[39:47])
            z = float(line[47:55])
            tmp = x,y,z    
            dt.append(tmp)
            a_n.append(sl[2])
        fname.close()
        return array(dt),a_n

    # Writes the transformed coordinates to file
    # Requires coordinates translated and rotated
    # Returns xyz file of rotated and translated
    def write_data(self,atomnames,coordinates):
        al = len(atomnames)
        #    assert al == len(coordinates),'Lengths of lists do not match'
        wname = open('superimposed.xyz','w')
        wname.write(str(al)+'\n Generated to Enzyme Design\n')
        for i in range(al):
            wname.write(atomnames[i]+'\t'+str(coordinates[i][0])+'   '+str(coordinates[i][1])+'   '+str(coordinates[i][2])+'\n')


    
    def write_pdb(self,atomnames,coordinates,fname='superimposed.pdb'):
        '''Write a pdb file of the aligned coordinates '''
        al = len(atomnames)
        # Name of chain
        ch = '  SUB'
        wname = open(fname,'w')
        for i in range(al):
            wname.write('ATOM   '+str(i).rjust(3)+'   '+atomnames[i]+' '+ch+' X '+str(i).rjust(3)+'    '+str(coordinates[i][0])[0:6]+'  '+str(coordinates[i][1])[0:6]+'  '+str(coordinates[i][2])[0:6]+'\n')
            
            

    # Superimpose a on b
    # Requires two arrays same length
    # Returns translation and rotation matrices
    def get_rotate_translate(self,a,b):
        assert shape(a) == shape(b), 'Matrix not of same dimensions'
        # Store number of rows
        nr = shape(a)[0]
        # Debug coordinates
        # print 'a coordinates', a
        # print 'b coordinates', b        
        # Get the normalized coordinates
        av_a = sum(a,axis=0)/nr
        av_b = sum(b,axis=0)/nr
        # Normalized matrices
        norm_a = a-av_a
        norm_b = b-av_b
        # Generate cross dispersion matrix for a and b
        cd = dot(transpose(norm_a),norm_b)
        # Do the SVD in order to get rotation matrix
        u,s,vt = linalg.svd(cd)
        # Rotation matrix
        rt_m = transpose(dot(transpose(vt),transpose(u)))
        # Due to reflection check determinant of matrix det A = -1
        if linalg.det(rt_m) < 0:
            vt[2] = -vt[2]
            rt_m = transpose(dot(transpose(vt),transpose(u)))
        tr_m = av_b-dot(av_a,rt_m)
        ## Debug
        ## print 'Determinant of alignment',linalg.det(rt_m)
        ##
        # Returning the translation and rotational matrices
        # print 'Rotation matrix',rt_m
        # print 'Translation matrix',tr_m
        return tr_m, rt_m

    # Returning the superimposed coordinates
    # Requires two arrays same length
    # Returns array of transformed coordinates
    def get_transformed_coor(self,a,b):
        tr, rt = self.get_rotate_translate(a,b)
        nw_coor = dot(a,rt) + tr
        return nw_coor


    # Computes the rmsd between aligned and original data set
    # Requires transformed array and original array
    # Returns rmsd 
    def get_rmsd(self,t_a,o_a):
        # Number of rows in matrices
        nr_r = shape(t_a)[0]
        # Difference between matrices
        df = (t_a - o_a) /nr_r
        return sqrt(sum(df*df))


    # Testing of methods
    def __test__():

        tt,an = get_data('coor.pdb')
        tmp = get_transformed_coor(tt,tt)
        print get_rmsd(tmp,tmp)
        
        print an
        
        write_data(an,tt)

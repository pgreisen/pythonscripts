from numpy import *

''' 

Superimpose one coordinate set on another
The algorithm is based on Landau mechanics
Returns the coordinates in pdb file

'''

class translate_rotate:

    # Superimpose a on b
    # Requires two arrays same length
    # Returns translation and rotation matrices
    def get_rotate_translate(self,a,b):

        assert shape(a) == shape(b), 'Array not of same dimensions'

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

        # Do the SVD to get rotation matrix
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

    # a contains the coordinates needed to be transformed

    def get_transformed_coor(self,a,b):
        tr, rt = self.get_rotate_translate(a,b)
        nw_coor = dot(a,rt) + tr
        return nw_coor


    def transform_ligand_coordinates(self,ligand_coordinates,translation_matrix,rotation_matrix):
        transformed_coordinates = {}
        for k,v in ligand_coordinates.items():
            transformed_coordinates[k] =  (dot(v,rotation_matrix) + translation_matrix )
        return transformed_coordinates
            

    # get rmsd between two dictionaries of atoms
    # Computes the rmsd between aligned and original data set
    # Requires transformed array and original array
    # Returns rmsd
    def get_dictionary_rmsd(self,native_coordinates,design_coordinates):
        number_of_atoms = 0
        df = 0.0
        for key in native_coordinates:
            if( key in design_coordinates ):
                number_of_atoms += 1
                x = native_coordinates[key] - design_coordinates[key]
                df += linalg.norm(x)
        return df/ number_of_atoms


    # Computes the rmsd between aligned and original data set
    # Requires transformed array and original array
    # Returns rmsd 
    def get_rmsd(self,t_a,o_a,t_atoms,o_atoms):
        # Number of rows in matrices
        nr_r = shape(t_a)[0]
        df = 0
        for i in range(nr_r):
            x = (t_a[i]-o_a[i])
            df = df + linalg.norm(x)
            # Debug
            # print round(df,3),t_atoms[i],o_atoms[i],round(linalg.norm(x),3)

        # Difference between matrices
        # df = (t_a - o_a) /nr_r
        return df/nr_r

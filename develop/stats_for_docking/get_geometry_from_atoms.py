#!/usr/bin/python
from math import sqrt
from numpy import *
import os,shutil,argparse
import pandas as pd
import shutil
'''
Fix_pdbfile get chain A as well as hetatms connected to chain A.
It takes all a pdb file as input and returns
the same file.

@parameter f: PDB file with metal site
@type       : string
@parameter m: Metal ion - default ZN 
@type       : string
@parameter n: Name of rosetta output pdb file 
@type       : string
@parameter l: Three ligand names with Zn and its first ligand
@type       : string
@parameter a: Coordinate file of aligned ligand and is append to the Rosetta input
@type       : string
@return     : constraint file, pdb file to use as input for rosetta


The default metal-protein ligand distance is set to 2.9 AA and for non-protein it is set to 2.7 AA

PTE_C23_PP_FINAL_DIMER_0001//A/ZNB`364/ZN

HO- - O;

'''


class GetGeometry:

    def __init__(self):
        # Parameters for the distance geometry interactions
        # Distance between metal and protein
        self.DISTANCEMETAL = 2.9
        # Distance between metal and hetero atom
        self.DISTANCEHET = 2.7

        self.metal_coor = {}

        self.METAL = "ZN"

        self.path = "./"

        self.coordination_sites_protein = {}

        self.chain = "A"
        self.geometry = {}
        self.metal_coordination = {}
        self.coordinates = {}
        self.name_of_pdb = ""
        # based on mespeus database with 3.0 Ang cutoff
        self.metal_distance_cutoff = 2.7
        self.restraint_cutoff = 3.0
        self.data_dir = ""


    def write_stats_file(self):
        with open("stats.dat", "w") as f:
            f.write("Metal site,residue,distance,angle,torsion\n")
            for key in self.metal_coordination:
                for key2 in self.metal_coordination[key]:
                    if(len(self.metal_coordination[key][key2]) == 3):
                        f.write(str(key)+","+str(key2)+","+self.metal_coordination[key][key2][0]+","+self.metal_coordination[key][key2][1]+","+self.metal_coordination[key][key2][2]+"\n")

                    else:
                        f.write(str(key)+","+str(key2)+","+self.metal_coordination[key][key2][0]+"\n")

    # Requires list
    # Returns one list
    def merge(self,seq):
        merged = []
        for s in seq:
            for x in s:
                merged.append(x)
        return merged

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


    # Requires list with pdb lines, metal coordinates ( set in the constructor )
    # Returns dictionary with ligands < DISTANCE from metal ion
    def get_protein_metal_geometry(self, PDB ):
        with open(PDB, 'r') as f:
            for line in f:
                if (line[21:22] == self.chain and line[0:4] == "HETA" and line[13:14] != 'H'):
                    tmp = line.split()
                    key = tmp[2] + "_" + tmp[3]
                    self.coordinates[key] = array([float(line[31:39]), float(line[39:47]), float(line[47:55])])

                if (line[21:22] == 'X' and line[0:4] == "HETA" and line[13:14] in ['O','P']):
                    tmp = line.split()
                    key = tmp[2] + "_" + tmp[3]
                    self.coordinates[key] = array([float(line[31:39]), float(line[39:47]), float(line[47:55])])

    def get_geometry(self,dictionary):
        for key in dictionary.keys():
            for key2 in dictionary.keys():
                if(key != key2):
                    dist = linalg.norm( dictionary[key]  -  dictionary[key2])
                    new_key = key+"_"+key2
                    self.geometry[self.name_of_pdb][new_key] = dist

    def write_stats_to_file(self):
        with open("stats.csv",'w') as f:
            f.write("PDB_name,dist_name,dist\n")
            for key in self.geometry.keys():
                for value in self.geometry[key].keys():
                    f.write(key+","+str(value)+","+str(round(self.geometry[key][value],2) )+"\n")

    def convert_to_csv(self):
        import pandas as pd
        df = pd.read_csv("stats.csv")
        df = df.pivot(index='PDB_name', columns='dist_name', values='dist')
        df.reset_index()
        df.to_csv("stats.csv")



    def analyse_geometry_pdbs(self):
        '''
        :return: 0
        '''
        pdbfiles = os.listdir('./')
        for pdbfile in pdbfiles:
            if (os.path.isfile(pdbfile) and str(pdbfile[-3:]) == 'pdb'):
                self.name_of_pdb = pdbfile.split('.')[0]
                self.geometry[self.name_of_pdb] = {}
                # set self.coordinates in constructor
                self.get_protein_metal_geometry(pdbfile)
                # set geometry dictionary in cnstr
                self.get_geometry(self.coordinates)
        self.write_stats_to_file()
        self.convert_to_csv()


    def merge_results(self,tagname):
        '''
        Merge the dock.sc file with the analysis of the pdb file
        :param tagname:
        :return:
        '''
        import glob
        scorefile =  glob.glob("*.sc")[0]
        df_dock = pd.read_csv(scorefile, skiprows=1, sep='\s+')
        df = pd.read_csv("stats.csv")
        new_df = pd.merge(left=df_dock, left_on='description', right=df, right_on='PDB_name')
        # filter docking poses based on distance restraint between oxygen phosphoryl and zinc B
        resname = tagname.split('_')[0]
        new_df = new_df[ (new_df['ZN_ZNB_O2_'+resname] <= self.metal_distance_cutoff ) | (new_df['ZN_ZNB_O1_'+resname] <= self.metal_distance_cutoff ) ]
        new_df = new_df[(new_df['atom_pair_constraint'] <= self.restraint_cutoff) & (new_df['coordinate_constraint'] <= self.restraint_cutoff)]
        new_df.to_csv(tagname+".csv")
        shutil.copy2(tagname+".csv", '../')


    def concatenate_df(self,resname):
        dfs = []
        drfiles = os.listdir('./')
        for csvfile in drfiles:
            if (os.path.isfile(csvfile) and str(csvfile[-3:]) == 'csv'):
                tmp_df = pd.read_csv(csvfile)
                dfs.append(tmp_df)
        total = pd.concat(dfs)
        total.to_csv(resname + "_total_data.csv")


    def collect_pdb_files_below_threshold(self,resname):
        import subprocess
        import glob
        list_of_pdbs = []
        df = pd.read_csv(resname + "_total_data.csv")
        # make directory

        cwd = os.getcwd()

        dst = cwd+"/"+resname+"_pdbs"
        os.makedirs(dst)
        # shutil.copy2(resname + "_total_data.csv", dst)

        for file in glob.glob(r'*.csv'):
            shutil.copy2(file, dst)

        for pdbfile in df['PDB_name']:
            tmpfile = pdbfile+".pdb"
            srcfile = self.find_file(tmpfile, './')
            self.copy_file(srcfile,dst )

        # tar zip directory and copy to dst
        self.make_tarfile(resname+"_pdbs.tgz" , dst)
        shutil.copy2(resname+"_pdbs.tgz", self.data_dir)


    def find_file(self,name, path):
        # find file in given path (the first match)
        for root, dirs, files in os.walk(path):
            if name in files:
                return os.path.join(root, name)


    def copy_file(self, src, dest):
        # copy file from source to destination
        try:
            shutil.copy(src, dest)
        except:
            print "err"
        #except shutil.Error as e:
        #    print 'Error: {0}' #.format(e)
        # eg. source or destination doesn't exist
        #except IOError as e:
        #print 'Error: {0}' # .format(e.strerror)

    def make_tarfile(self, output_filename, source_dir):
        import tarfile, os
        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

    def main(self):
        parser = argparse.ArgumentParser(\
            description="Compute the Pareto front between the list of data under construction ")
        parser.add_argument("--ylabel", dest="ylabel", help="ylabel", default="IFE")

        input_variables = parser.parse_args()
        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # get name of working directory
        workdir = str(os.getcwd()).split('/')[-1]

        # MKDIR
        self.data_dir = os.getcwd()+"/"+workdir+"_data"
        os.makedirs(self.data_dir)

        # structure of directories SVR etc.
        # and subdictory contains 1-10
        drcts = os.listdir('./')
        for dr in drcts:
            if (os.path.isdir(dr) and len(dr) == 3):
                # changing to subdirectory
                os.chdir(dr)
                print "Working on this directory: ", dr
                subdrs = os.listdir('./')
                for subdr in subdrs:
                    if(os.path.isdir(subdr) and subdr != "init"):
                        os.chdir(subdr)
                        # generate csv and merge with docking output
                        self.analyse_geometry_pdbs()
                        tagname = str(dr)+"_"+subdr
                        # merge results in each directory between pdb and dock output
                        self.merge_results(tagname)
                        os.chdir('../')

                # reset value after each ligand
                self.geometry.clear()
                self.coordinates.clear()

                self.concatenate_df(str(dr))
                # get relevant pdb files
                self.collect_pdb_files_below_threshold( str(dr) )

                os.chdir("../")


if __name__ == "__main__":
    run = GetGeometry()
    run.main()

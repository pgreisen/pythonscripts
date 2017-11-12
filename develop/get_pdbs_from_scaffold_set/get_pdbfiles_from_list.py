'''

Getting files from the /lab/shared/scaffolds/
as there are different ligands in the set the user has to specify 
which ligand to get e.g.,

number = '1'
pos_zero = '1.pos'
gridlig_zero = '1.gridlig'

1. The files are copied from
   /lab/shared/scaffolds/ with the $basename/$pdbname/

2. Take a lists of pdbs and copies each pdb into a directory with its
   position file present for the biggest ligand. The lists have to be constructed with the
   pdb name as the first 4 characters.

   python get_pdblist.py list_of_pdbs

'''

import sys, shutil, os, subprocess

list_of_pdbs = sys.argv[1]
path_relax = '/lab/shared/scaffolds/'


def read_inputlist(list_of_pdbs):
    pdbfl = open(list_of_pdbs,'r')
    pdbfiles = []
    for line in pdbfl:
        pdbfiles.append(str(line[0:4]))
    return pdbfiles

def get_gridlig_pos_path(pdb_path):
    print pdb_path

# string of position ( 0.pos or 1.pos )
# string of gridlig
def gridlig_present_is_present(pdb_path,position,gridlig):

    is_there_pos = False
    is_there_gridlig = False

    files = os.listdir(pdb_path)

    for fl in files:

        if(str(fl[-5:]) == position):
            is_there_pos = True

        if(str(fl[-9:]) == gridlig):
            is_there_gridlig = True

    return is_there_pos, is_there_gridlig


def generate_dir_and_copy_files(arg1,arg2,arg3,directoryname):

    try:

        subprocess.Popen(arg1,shell=True)
        subprocess.Popen(arg2,shell=True)
        subprocess.Popen(arg3,shell=True)
        shutil.os.mkdir(directoryname)

    except:
        print "Directory was not made and files not copied"

def main():
    l_o_pdbs = read_inputlist(list_of_pdbs)
    number = '1'
    pos_zero = '1.pos'
    gridlig_zero = '1.gridlig'

    for elem in l_o_pdbs:
        pdb_path = str(path_relax)+str(elem[1:3])+'/'+str(elem)+'/'

        if(os.path.exists(pdb_path)):

            grid,pos = gridlig_present_is_present(pdb_path,pos_zero,gridlig_zero)

            if(grid == True and pos == True):
                
                gridlig_file = pdb_path+'/'+elem+'_nohet_1.pdb_'+gridlig_zero
                pos_file = pdb_path+'/'+elem+'_nohet_1.pdb_'+pos_zero
                pdb_file = pdb_path+'/'+elem+'_nohet_1_relax.pdb'


                if(os.path.exists(gridlig_file) and os.path.exists(pos_file) and os.path.exists(pdb_file)):
                    dst = elem+'_'+number
                    arg_grid = 'cp '+pdb_path+'*nohet_1.pdb_'+gridlig_zero+' ./'+dst+'/ligand.gridlig'
                    arg_pos = 'cp '+pdb_path+'*nohet_1.pdb_'+pos_zero+' ./'+dst+'/protein.pos'
                    arg_pdb = 'cp '+pdb_path+'/*_nohet_1_relax.pdb ./'+dst+'/'

                    generate_dir_and_copy_files(arg_grid,arg_pos,arg_pdb,dst)

                else:
                    print 'one of the files did not exist (maybe the 0 file?) ',gridlig_file
        else:
            print "This path does not exist",pdb_path

if __name__ == "__main__":
        main()

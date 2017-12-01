from GetLigandPositions import *
import argparse,shutil,os,subprocess,re

'''
It takes a pdbfile with a ligand present and generates a lig-file
for matchers. 

@requires: pdbfile, position file, and ligand name


'''

PATH='/work/greisen/rosetta/develop_rosetta/bin/'
EXE=''

# @require pdbfile list
# @return list with ligand or protein
def get_xyz(pdbfile,atomtype):

    tmp_list = []    
    for line in pdbfile:
        if line[0:4] == atomtype:
            tmp_list.append(line)
    return tmp_list


# @require list, name
# @write file with name
def write_outfile(namelist,outname):
    tmp_write = open(outname,'w')
    for line in namelist:
        tmp_write.write(line)
    tmp_write.close()

# @require to be in a directory with the
# ligand and protein files
# @generates the gridlig file
def run_grid(tmppr,tmplg):
    exe = '~/rosetta/rosetta_source/bin/gen_lig_grids.linuxiccrelease -s '+tmppr+' '+tmplg+' -database ~/rosetta/rosetta_database @../../InputFiles/gen_lig.flags'
    subprocess.Popen(exe,shell=True).wait()


# @require pdbfile and the atom names
# @generates a position file with distances less 6 Aangstrom
def gen_pos(pdbfile,atomnames):
    exe = 'python ~/pythonbin/getligandpositionmatchermain.py -f '+pdbfile+' -n '+atomnames
    subprocess.Popen(exe,shell=True).wait()


# @requires: pdbfile, position file, and ligand name
def setup(pdbfile):

    gp = GetLigandPositions()

    # Get pdbfile
    p_pdbfile = gp.get_pdbfile(pdbfile)

    protein = get_xyz(p_pdbfile,'ATOM')
    ligand = get_xyz(p_pdbfile,'HETA')

    # Variables used
    tmppr = 'tmpprotein.pdb'
    tmplg = 'tmpligand.pdb'

    write_outfile(protein,tmppr)
    write_outfile(ligand,tmplg)

    run_grid(tmppr,tmplg)

    # input for aromatic system where a virtual atom has not
    # been specified
    gen_pos(pdbfile,'N1 O1 aro-C8-C19 aro-C10-C22')


# removes all the temporary files such as protein.pdb, ligand.pdb, and protein.pdb_0.pos
def clean():
    exe = 'rm -f tmpprotein.pdb tmpligand.pdb tmpprotein.pdb_0.pos'
    subprocess.Popen(exe,shell=True).wait()


def main():
    path = './'

    drs = os.listdir(path)
    for dr in drs:
        if os.path.isdir(dr):
            os.chdir(dr)

            files = os.listdir(path)


            for fl in files:
                if os.path.isfile(fl):

                    evl = re.match('(.*)_00[0-9]*.pdb',fl)

                    
                    if(evl != None):
                        setup(fl)

            clean()
            os.chdir('../')
            print 'Done with directory ',dr



if __name__ == "__main__":
    main()    
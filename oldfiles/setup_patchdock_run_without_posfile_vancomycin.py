from optparse import OptionParser
import os,shutil,subprocess

'''

Runs patchdock on a set of pdbs in a directory, which only contains ONE pdbfile and ONE
position file.

python setup_patchdock_run.py -l ligandpdb -s scaffoldset

All the scaffold set are located in /work/greisen/Scaffolds/


# needs to be added - status report have all the job finished correctly
# multi-threading of jobs. 

'''

def make_output_dir(ligandname,ligandpdb):
    os.mkdir(ligandname)
    shutil.move(ligandpdb,ligandname)


def copy_scaffold_dirs(scaffoldpath,abs_path):
    dirs = os.listdir(scaffoldpath)
    for dir in dirs:
        if os.path.isdir(scaffoldpath+'/'+dir):
            shutil.copytree(scaffoldpath+'/'+dir,abs_path+'/'+dir)

def run_patchdock(ligandname,pdbfile,positionfile):
    debug_executable = open('debug_exe.sh','w')
    exe = '/work/greisen/files/PatchdockExecutable/make_patchdock_vancomycin.pl'

    exe_run = 'perl '+exe+' '+ligandname+' '+pdbfile
    debug_executable.write(exe_run)
    subprocess.Popen(exe_run,shell=True).wait()
    debug_executable.close()

def setup_run(ligand_directory,ligandname):
    pdbfile = ''
    positionfile = ''
    path = './'
    dirs = os.listdir(path)
    for dir in dirs:
        if os.path.isdir(dir) and dir !=ligandname:
            os.chdir(dir)
            files = os.listdir(path)
            for fl in files:
                print fl
                if fl.endswith('pdb'):
                    # The patchdock executable is hardcoded with this
                    # should be changed at some point
                    pdbfile = fl.split('.')[0]

                elif fl.endswith('pos'):
                    positionfile = fl

                else:
                    continue
            if pdbfile is not None:
                run_patchdock(ligand_directory,pdbfile,positionfile)
            os.chdir('../')




def main():

    parser = OptionParser()

    parser.add_option('-l',dest='ligandpdb',help='Ligand pdb file')
    parser.add_option('-s',dest='scaffoldset',default='Moad',help='Scaffold set used for the patchdocking')
    parser.add_option('-n',dest='ligandname',default='ligand1',help='Name of the ligand')

    (options,args) = parser.parse_args()

    ligandpdb = options.ligandpdb

    ligandname = ligandpdb.split('.')[0]

    scaffoldname = options.scaffoldset

    scaffoldpath = '/work/greisen/Scaffolds/'+scaffoldname

    # Ligand directory absolute path of ligand pdb
    # without extension
    abs_path = os.path.abspath('./')
    ligand_directory = str(abs_path)+'/'+ligandname+'/'+ligandname


    # Make output directory with ligand file
    make_output_dir(ligandname,ligandpdb)

    # Copy directories to working directory
    copy_scaffold_dirs(scaffoldpath,abs_path)

    # Sets up directories and run them
    setup_run(ligand_directory,ligandname)


    
if __name__ == '__main__':
    main()

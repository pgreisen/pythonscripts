import sys, shutil, os, subprocess
'''

Setup parallel run on the digs

'''

def get_native_pdbfile( pdbfile ):
    pdbname = pdbfile.lower()

    pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+pdbname+'_nohet_1_relax.pdb'
    # move_files = 'scp dig1:'+pdb_path+' .'
    # 
    # subprocess.Popen(move_files,shell=True).wait()
    return pdb_path


print "1: run string e.g ~/files/SecondaryOrientationFilter/run.sh "
print "2: parameter file ( full path )"
print "3: Number of processors"
print "4: the number of the split-string to keep for identifying the pdb"


run_string = sys.argv[1]
number_of_processors = sys.argv[2]


path = './'
dirs = os.listdir( path )

run_file = open( "run_file.txt", 'w' )

for dir in dirs:
    if ( os.path.isdir ( dir ) ):
        os.chdir( dir )
        pdbs = os.listdir(".")

        run_file.write( "sh "+str(run_string)+" \n"  )
        os.chdir("../")



exe = "cat run_file.txt | parallel -j "+str( number_of_processors )
print exe
subprocess.Popen(exe, shell=True).wait()

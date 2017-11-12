import sys, shutil, os, subprocess
'''
Setup parallel run on the digs

'''

def get_native_pdbfile( pdbfile ):
    pdbname = pdbname.lower()

    pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+pdbname+'_nohet_1_relax.pdb'

    # move_files = 'scp dig1:'+pdb_path+' .'
    # 
    # subprocess.Popen(move_files,shell=True).wait()
    return pdb_path


print "1: run string e.g ~/files/SecondaryOrientationFilter/run.sh "
print "2: parameter file ( full path )"
print "3: Number of processors"

run_string = sys.argv[1]
parameterfile = sys.argv[2]
number_of_processors = sys.argv[3]

path = './'
dirs = os.listdir( path )

run_file = open( "run_file.txt", 'w' )

for dir in dirs:
    if ( os.path.isdir ( dir ) ):
        run_file.write( "sh "+str(run_string)+" "+str(dir)+"/*.pdb "+str( parameterfile )+"\n"  )

exe = "cat run_file.txt | parallel -j "+str( number_of_processors )
print exe
subprocess.Popen(exe, shell=True).wait()

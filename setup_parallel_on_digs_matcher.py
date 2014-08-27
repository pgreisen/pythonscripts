import sys, shutil, os, subprocess
'''
Setup parallel run on the digs for matching

'''

print "1: run string e.g ~/files/SecondaryOrientationFilter/run.sh "
print "2: parameter file ( full path )"

print "3: Name of residue e.g. FEN"
print "4: Cstfile - hardcoded is the location of it"
print "5: Number of processors"

run_string = sys.argv[1]
parameterfile = sys.argv[2]
residue = sys.argv[3]
cstfile = sys.argv[4]
number_of_processors = sys.argv[5]


path = './'
dirs = os.listdir( path )

run_file = open( "run_file.txt", 'w' )

for dir in dirs:
    if ( os.path.isdir ( dir ) ):
        run_file.write( "cd "+dir+"/; sh "+str(run_string)+" *.pdb "+str( parameterfile )+" "+str(residue)+"/work/greisen/Projects/LigandBinding/VX/Cstfiles/"+str(cstfile)+" > logfile_"+str(cstfile)+";  cd ..;\n"  )

exe = "cat run_file.txt | parallel -j "+str( number_of_processors )
print exe
subprocess.Popen(exe, shell=True).wait()

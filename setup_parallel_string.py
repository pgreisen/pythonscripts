import sys, shutil, os, subprocess
'''
Setup parallel run on the digs

'''
print "1: run string e.g ~/files/SecondaryOrientationFilter/run.sh "
print "2: Number of processors"

run_string = sys.argv[1]
number_of_processors = sys.argv[2]

path = './'
dirs = os.listdir( path )

outfile = "run_file_sty_1_2.txt"

run_file = open( outfile, 'w' )

for dir in dirs:
    if ( os.path.isdir ( dir ) ):
        run_file.write( "cd "+dir+"/; "+str(run_string)+";   cd ..;\n"  )

exe = "cat "+str(outfile)+" | parallel -j "+str( number_of_processors )
print exe
subprocess.Popen(exe, shell=True).wait()

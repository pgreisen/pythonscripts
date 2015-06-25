import sys

pdbfile = sys.argv[1]

pdbname = sys.argv[2]

newfile = open(pdbname,'w')

dummy = 1

with open(pdbfile,'r') as fl:
    for line in fl:
        if( len(line.strip()) == 67 ) :
            newline = line[0:6]+line[8:]
            newfile.write( newline )
        else:
            newfile.write( line )

        if(dummy % 17 == 0):
            newfile.write("END\n")
            #print dummy
        dummy = dummy +  1


import os,sys

startposition = int( sys.argv[1] )
endposition = int( sys.argv[2] )

tmpfile = open("pos.pos",'w')

while startposition <= endposition:
    tmpfile.write(str(startposition)+" ")
    startposition += 1

tmpfile.close()

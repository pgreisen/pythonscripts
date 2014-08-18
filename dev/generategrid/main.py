import sys
from numpy import *

def get_coordinates(inputfile):

    coordinates = []

    with open(inputfile) as f:
        for line in f:
            x = str(line[30:38]).rstrip()
            y = str(line[38:46]).rstrip()
            z = str(line[46:54]).rstrip()
            coordinates.append( array( [ float(x), float(y), float(z) ] ) )
    print "coordinates are ", coordinates
    print "mean of the coordinates are ", mean(coordinates, axis=0)
    return coordinates


def generate_grid_files( coordinates ):


    tmpfile = open("gridfile.dat", 'w')

    gridpoints = linspace(-12,12,64)

    geometrical_center =  mean(coordinates, axis=0)

    slice = 1

    for x in  gridpoints:
        tmpfile.write("\n Grid "+str(slice)+"\n")
        slice = slice + 1
        for y in gridpoints:
            tmpfile.write('\n')
            for z in gridpoints:
                tmp_vector = array([x,y,z])+geometrical_center
                on_grid = False
                for i in coordinates:
                    if( linalg.norm( i - tmp_vector ) < 1.4 ):
                        on_grid = True
                        continue
                if( on_grid ):
                    tmpfile.write( "1 " )
                else:
                    tmpfile.write( "0 " )
    print "done"

    tmpfile.close()

    #gridsize = 0.375
    #gridpoints = 64

    #with open("gridfile.dat", 'w') as f:
    #
    #
    #    if ():
    #
    #    else:
    #        f.write( 1 )



def main():
    inputfile = sys.argv[1]
    coordinates = get_coordinates( inputfile )

    generate_grid_files( coordinates )


if __name__ == "__main__":
    main()


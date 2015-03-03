import sys, shutil, os, subprocess, argparse


def get_file(cstfile):
    indexes = []
    with open(cstfile,'r') as f:
        for line in f:
            if(line[0] == '1'):
                indexes.append( line.split()[1] )
                indexes.append( int(line.split()[1]) -1  )
                indexes.append( int( line.split()[1] ) +1 )

    return indexes


def write_index_to_file(index):
    with open("position_patchdock",'w') as f:
        for ndx in index:
            f.write(str(ndx) +" "+"A\n")

def main():
    cstfile = sys.argv[1]

    index = get_file(cstfile)
    write_index_to_file(index)


if __name__ == "__main__":
   main()



#
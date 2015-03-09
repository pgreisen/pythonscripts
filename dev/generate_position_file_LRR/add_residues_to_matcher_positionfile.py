import sys, shutil, os, subprocess, argparse


def get_file(cstfile):
    indexes = []
    with open(cstfile,'r') as f:
        for line in f:
            indexes.append( line )
    return indexes


def write_index_to_file(index,additionalpositions,nm):
    dummy = 0
    with open( str(nm)+".pos",'w') as f:
        for ndx in index:
            if( dummy < 2 ):
                f.write(str(ndx) )
                dummy = dummy +  1
        f.write("2: "+additionalpositions+"\n")
        f.write("3: "+additionalpositions+"\n")
        f.write("4: "+additionalpositions+"\n")




def get_nearby_residues(additional_file):

    with open(additional_file,'r') as f:
        for line in f:
            add_positions = line
    return add_positions


def main():
    posfile = sys.argv[1]
    nm = sys.argv[2]
    additional_file = sys.argv[3]

    add_residues = get_nearby_residues(additional_file)

    index = get_file(posfile)
    write_index_to_file(index, add_residues,nm)


if __name__ == "__main__":
   main()



#
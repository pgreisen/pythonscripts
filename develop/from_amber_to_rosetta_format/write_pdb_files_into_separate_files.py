from pdbfile import *




def main():
    pfile = pdbfile()
    file = sys.argv[1]
    pdbfile = pfile.read_file(file)
    pfile.split_into_multiple_pdb_files(pdbfile,'MODEL')
    

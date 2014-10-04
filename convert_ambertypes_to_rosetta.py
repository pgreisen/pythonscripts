import sys

def convert_amber_atomtype_to_rosetta_atomtype(amberfile):
    """
    @requieres that ligand_am1_bcc.mol2 files is already generated
    
    @returns a mol2 with atomtypes readable for Rosetta's mol_to_parameter scripts
    
    """

    tmpfile = open("tmp.mol2", 'w')
    with open(amberfile,'r') as f:
        atoms = False
        for line in f:
            
            if ( len(line) > 13 and line.find("@<TRIPOS>ATOM") >-1.0):
                atoms = True
            elif ( len(line) > 13 and line.find("@<TRIPOS>BOND") >-1.0):
                atoms = False
                
            elif( atoms == True and len(line) > 75 ):
                tmp_characters = line[47]+"."+line[48]
                line = line[0:47]+tmp_characters+line[50:]
            tmpfile.write(line)
    tmpfile.close()

def main():
    amberfile = sys.argv[1]
    convert_amber_atomtype_to_rosetta_atomtype(amberfile)

if __name__ == "__main__":
   main()

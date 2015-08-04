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
                atms =line[47:49].upper()

                if(atms[0] == 'H' ):
                    tmp_characters = line[46]+atms[0]+" "+line[49]
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'CC' ):
                    tmp_characters = line[46]+"C.2"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'CE' ):
                    tmp_characters = line[46]+"C.2"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'CF' ):
                    tmp_characters = line[46]+"C.2"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'SS' ):
                    tmp_characters = line[46]+"S.3"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'OS' ):
                    tmp_characters = line[46]+"O.2"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'O ' ):
                    tmp_characters = line[46]+"O.3"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'OH' ):
                    tmp_characters = line[46]+"O.3"
                    line = line[0:46]+tmp_characters+line[50:]

                elif(atms == 'p5' ):
                    tmp_characters = line[46]+"P.3"
                    line = line[0:46]+tmp_characters+line[50:]

                else:
                    tmp_characters = line[46]+atms[0]+"."+line[48]
                    line = line[0:46]+tmp_characters+line[50:]

            tmpfile.write(line)
    tmpfile.close()

def main():
    amberfile = sys.argv[1]
    convert_amber_atomtype_to_rosetta_atomtype(amberfile)

if __name__ == "__main__":
   main()

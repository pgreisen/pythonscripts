import sys


class convert_md_results_into_rosetta_types:

    def readfile(self,file):
        tmp_file = open(file,'r')
        pdbfile = tmp_file.readlines()
        tmp_file.close()
        assert isinstance(pdbfile, object)
        return pdbfile


    # replace residue name in pdb file format
    # and exchange ATOM  to HETATM
    def replace_residue_name(self,new_residue_name,old_residue_name,pdbfile):
        new_pdbfile = []
        for line in pdbfile:
            if line[17:20] == old_residue_name:
                line = 'HETATM'+line[6:17]+new_residue_name+line[20:]
                new_pdbfile.append(line)


            elif line[17:20] == "ARN":
                line = line[0:17]+"ARG"+line[20:]
                new_pdbfile.append(line)



            elif line[17:20] == "LYN":
                line = line[0:17]+"LYS"+line[20:]
                new_pdbfile.append(line)


            elif line[17:20] == "ASH":
                line = line[0:17]+"ASP"+line[20:]
                new_pdbfile.append(line)

            elif line[17:20] == "HOH":
                print line
                continue


            elif line[17:20] == "GLH":
                line = line[0:17]+"GLU"+line[20:]
                new_pdbfile.append(line)

            else:
                new_pdbfile.append(line)


        return new_pdbfile


    # write file to disk
    def write_file(self,new_pdbfile):
        new_file = open('ligand_from_md.pdb','w')
        for line in new_pdbfile:
            new_file.write(line)


    def main(self):
        print "The script is executed like this:\n 1. Inputfile \n 2. New residue name \n 3. Old residue name\n"
        inputfile = sys.argv[1]
        new_residue_name = sys.argv[2]
        old_residue_name = sys.argv[3]

        pdbfile = self.readfile( inputfile )
        new_pdbfile = self.replace_residue_name(new_residue_name,old_residue_name,pdbfile)
        self.write_file(new_pdbfile)
        print "Done"



if __name__ == "__main__":
   run = convert_md_results_into_rosetta_types()
   run.main()


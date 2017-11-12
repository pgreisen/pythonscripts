class SetupSolventMinimization():

    def write_parameter_file(self,residuenumber):
        template = '''
initial minimization solvent + ions
&cntrl
imin   = 1,
maxcyc = 20000,
ncyc   = 10000,
ntb    = 1,
ntr    = 1,
cut    = 10.0
/
Hold the protein and ligand fixed
500.0
RES 1 '''+str(residuenumber)+'''
END
END
    '''
        return template

    def get_number_of_residues(self,pdbfile):
        tmpfile = open(pdbfile,'r')
        max_number = 0
        for line in tmpfile:
            try:
                query_number = int( line[23:26].strip() )

                if( query_number > max_number):
                    max_number = int(line[23:26])
            except:
                pass
        return max_number

    def write_out_parameterfile(self, list):
        with open("minimization_of_solvent.in",'w') as f:
            for line in list:
                f.write(line)


    def main(self,pdbfile):
        
        max_number_residues = get_number_of_residues(pdbfile)

        list_template = write_parameter_file(max_number_residues)

        write_out_parameterfile(list_template)


#if __name__ == "__main__":
#    setup = SetupSolventMinimization()
#    setup.main(pdbfile)

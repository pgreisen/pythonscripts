'''

This is a helper function to generate the script for analysing the rotamers from the MD simulation.

amber is executed like this

~/amber12/bin/ptraj SRO_FILE.prmtop < rotamer_analysis.ptraj


'''

class AmberAnalysisFormat:

    def __init__(self):
        pass


    def number_of_trajecties(self,number_of_directories):
        template_string = ""
        assert isinstance(number_of_directories, int)
        i = 1
        while i <= number_of_directories:

            template_string = template_string+"trajin MD_"+str(i)+"/mdcrd\n"
            i += 1
        return template_string


    def get_chi_angles(self,key_to_dictionary):
        """

        @rtype : list
        """
        residue_atom_names = {
        'ARG' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','NE'],['CG','CD','NE','CZ'],['CD','NE','CZ','NH1'] ],
        'ASN' : [ ['N','CA','CB','CG'],['CA','CB','CG','OD1']],
        'ASP' : [ ['N','CA','CB','CG'],['CA','CB','CG','OD1']],
        'CYS' : [ ['N','CA','CB','SG'] ],
        'GLN' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','OE1']],
        'GLU' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','OE1']],
        'HIS' : [ ['N','CA','CB','CG'],['CA','CB','CG','ND1'] ],
        'ILE' : [ ['N','CA','CB','CG1'],['CA','CB','CG1','CD1']],
        'LEU' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD1'] ],
        'LYS' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','CE'],['CG','CD','CE','NZ']],
        'MET' : [ ['N','CA','CB','CG'],['CA','CB','CG','SD'],['CB','CG','SD','CE'] ],
        'PHE' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD1'] ],
        'PRO' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD'], ['CB','CG','CD','N'],['CG','CD','N','CA'] ],
        'SER' : [ ['N','CA','CB','OG']],
        'THR' : [ ['N','CA','CB','OG1'] ],
        'TRP' : [ ['N','CA','CB','CG'], ['CA','CB','CG','CD1'] ],
        'TYR' : [ ['N','CA','CB','CG'],['CA','CB','CG','CD1'] ],
        'VAL' : [ ['N','CA','CB','CG1'] ]
        }
        return residue_atom_names[ key_to_dictionary ]

    def amber_template(self, nof, number_of_residues, dihedral_string ):
        template = '''
# load the trajectory
# number of trajectories
# loop over number of trajectories
'''+nof+'''

# Remove water
strip :WAT

# need rms to remove rotation and translation
rms first mass :1-'''+str(number_of_residues)+'''
# calculate the dihedral angles for the main interacting residues
'''+dihedral_string+'''

        '''
        #print template
        return template




    def get_number_of_residues_pdb(self,pdbfile):
        i = 1
        for line in pdbfile:
            if( line[0:4] == "ATOM" ):
                tmp_residue = int(line[22:26].strip())
                if( tmp_residue > i ):
                    i = tmp_residue
        return i



    def get_dihedral_string(self,protein_positions):

        template_string = ''

        dummy = 1

        for key,value in protein_positions.iteritems():

            tmp_list = self.get_chi_angles ( protein_positions[key] )

            for j in tmp_list:

                tmp = "dihedral "+str( value )+"_"+str( key )+" :"+str(key)+"@"+j[0]+" :"+str(key)+"@"+j[1]+" :"+str(key)+"@"+j[2]+" :"+str(key)+"@"+j[3]+" out "+str( value )+"_"+str( key )+"_"+str(dummy)+".dat \n"

                template_string = template_string+tmp

                dummy += 1

        # print template_string

        return template_string


    def write_amber_analysis_to_file(self,protein_positions,pdbfile):

        dihedral_string = self.get_dihedral_string( protein_positions )

        number_of_directories = 200

        nof = self.number_of_trajecties(number_of_directories)

        number_of_residues = self.get_number_of_residues_pdb( pdbfile )

        am_file = self.amber_template( nof, number_of_residues, dihedral_string )


    def main(self):
        print "hello"



if __name__ == "__main__":
   run = AmberAnalysisFormat()
   run.main()

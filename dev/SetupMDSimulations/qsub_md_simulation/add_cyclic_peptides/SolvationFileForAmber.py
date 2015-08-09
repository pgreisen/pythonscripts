import os, shutil, sys, subprocess
# from numpy import *

class SolvationFileForAmber:

    def get_template(self,disulfide_pairs):
        ds_string = ""
        for i in disulfide_pairs:
            pair = i.split(',')

            ds_string = ds_string +"bond PRT."+str(pair[0])+".SG PRT."+str(pair[1]).strip()+".SG\n"
        
        
        template = '''
source leaprc.gaff

#loadoff phosphoaa10.lib
loadoff ions08.lib

loadamberparams frcmod.ionsjc_tip3p

PRT = loadpdb minimized.pdb

solvateBox PRT TIP3PBOX 10
'''+ds_string+'''
addions PRT Na+ 0
addions PRT Cl- 0

saveamberparm PRT solv.prmtop solv.inpcrd

quit
'''
        with open("parameters_for_solvation.sh",'w') as f:
            f.write(template)


        return template

    def get_template_ligand(self,ligandname,parm,lib):
        template = '''
source leaprc.gaff
loadoff ions08.lib
loadoff '''+lib+'''
loadamberparams frcmod.ionsjc_tip3p
loadamberparams '''+parm+'''
PRT = loadpdb minimized.pdb
solvateBox PRT TIP3PBOX 10
addions PRT Na+ 0
addions PRT Cl- 0
saveamberparm PRT solv.prmtop solv.inpcrd
quit'''
        with open("parameters_for_solvation.sh",'w') as f:
            f.write(template)

    def get_template_protein_ligand_complex(self,ligandname,parm,lib,disulfide_pairs):
        ds_string = ""
        for i in disulfide_pairs:
            pair = i.split(',')
            ds_string = ds_string +"bond PRT."+str(pair[0])+".SG PRT."+str(pair[1]).strip()+"SG\n"
        template = '''
source leaprc.gaff
loadoff ions08.lib
loadoff '''+lib+'''
loadamberparams frcmod.ionsjc_tip3p
loadamberparams '''+parm+'''
PRT = loadpdb minimized.pdb
'''+ds_string+'''
solvateBox PRT TIP3PBOX 10
addions PRT Na+ 0
addions PRT Cl- 0

saveamberparm PRT solv.prmtop solv.inpcrd

quit'''
        with open("parameters_for_solvation.sh",'w') as f:
            f.write(template)


    def get_disulfide_pairs(self):
        tmpfile = open('../Minimization/disulfide_pairs.txt','r')
        disulfide_paris = []
        for i in tmpfile:
            disulfide_paris.append(i)
        return disulfide_paris


    def get_cyclic_pairs(self):
        tmpfile = open("../Minimization/first_last.txt",'r')
        cyclic_paris = []
        for i in tmpfile:
            cyclic_paris.append(i)
        return cyclic_paris


    # def get_template_cp(self,disulfide_pairs):
    def get_template_cp(self, disulfides, ligname,parmfile,libfile ):
        ds_string = ""
        for i in disulfide_pairs:
            pair = i.split(',')

            ds_string = ds_string +"bond PRT."+str(pair[0])+".SG PRT."+str(pair[1]).strip()+".SG\n"

        cp = self.get_cyclic_pairs()

        bond = "bond PRT."+str( str( cp[0]).strip() )+".N PRT."+str( str( cp[1]).strip() )+".C\n"

        template = '''
source leaprc.gaff
#loadoff phosphoaa10.lib
loadoff ions08.lib
loadamberparams frcmod.ionsjc_tip3p
PRT = loadpdb minimized.pdb
loadoff ions08.lib
loadoff '''+lib+'''
loadamberparams frcmod.ionsjc_tip3p
loadamberparams '''+parm+'''

# Add cyclic bond
'''+bond+'''

solvateBox PRT TIP3PBOX 10
'''+ds_string+'''
addions PRT Na+ 0
addions PRT Cl- 0

saveamberparm PRT solv.prmtop solv.inpcrd

quit
'''
        with open("parameters_for_solvation.sh",'w') as f:
            f.write(template)

        return template



    def main(self,cyclic_peptide=False):

        disulfide_pairs = self.get_disulfide_pairs()

        if( cycle_peptide ):
            amber_solvation = self.get_template_cp(disulfide_pairs)

        else:                        
            amber_solvation = self.get_template(disulfide_pairs)

        ds_file = open("parameters_for_solvation.sh",'w')

        for line in amber_solvation:
            ds_file.write(line)


if __name__ == "__main__":
    ds = SolvationFileForAmber()
    ds.main()

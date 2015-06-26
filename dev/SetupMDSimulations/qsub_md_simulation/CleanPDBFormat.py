import os, shutil, sys, subprocess
from numpy import *

'''

@requires:
        libraryfile, parameterfile, pdbfile
@returns:
        clean pdbfile, histidine protonation states set according to amber,
        insert disulfides, generates amber parameter file
        
'''

class CleanPDBFormat:

    def __init__(self):
        self.DISULFIDESDISTANCE = 2.5


    def get_atom_numbers_for_cyclic_peptide_bond(self):
        '''
        


        '''
        pdbfile = open("design1.pdb", 'r')
        first_atom = ""
        last_atom = ""
        for line in pdbfile:
            residue_nr = str(line[22:26]).strip()
            atom_name = str(line[13:15]).strip()

            if( residue_nr == self.first_residue and atom_name == 'N' ):
                first_atom = str(line[7:11]).strip()
                # print "FIRST ATOM DEBUG",line

            elif( residue_nr == self.last_residue and atom_name == 'C' ):
                last_atom = str(line[7:11]).strip()
                # print "SECOND ATOM DEBUG",line
        return first_atom, last_atom



    def set_correct_histidines(self,tmplist, hie_numbers):
        '''

        '''
        correct_list = []

        for line in tmplist:

            if (line[17:19] == "HI" and line[23:26].strip() in hie_numbers):

                line = str(line[0:17]) + "HIE" + str(line[20:])

                correct_list.append(line)

            elif (line[17:19] == "HI"):

                line = str(line[0:17]) + "HID" + str(line[20:])

                correct_list.append(line)


            else:
                correct_list.append(line)

        return correct_list


    def set_disulfide_pairs(self,tmplist, residuepairs):
        # mega ugly hack
        residuenumbers = []
        for i in residuepairs:
            tmp = i.split(',')
            residuenumbers.append(tmp[0])
            residuenumbers.append(tmp[1])
        
        correct_list = []
        for line in tmplist:

            if (line[17:19] == "CY" and line[23:26].strip() in residuenumbers):

                line = str(line[0:17]) + "CYX" + str(line[20:])
                correct_list.append(line)

            else:
                correct_list.append(line)

        return correct_list


    def get_sulfur_atom_coordinates(self,pdbfile):
        sulfur_coordinates = {}

        for line in pdbfile:

            if ( line[0:4] == 'ATOM' and line[13:14] == 'S'):
                sulfur_coordinates[line[22:26].strip()] = array(
                [float(line[31:39]), float(line[39:47]), float(line[47:55])])

        return sulfur_coordinates


    def get_disulfide_pairs(self,sulfur_coordinates):
        """
        :param sulfur_coordinates:
        :return: pairs of disulfides
        """
        disulfides_pairs = []

        tmpfile = open("disulfide_pairs.txt",'w')

        already_done = []

        for k,v in sulfur_coordinates.iteritems():
            for kk,vv in sulfur_coordinates.iteritems():

                if ( kk != k) and kk not in already_done:
                    if linalg.norm( v - vv ) <= self.DISULFIDESDISTANCE:
                        tmpfile.write(k+','+kk+"\n")
                        already_done.append(k)
                        pair = k+','+kk
                        disulfides_pairs.append(pair)
                    else:
                        continue
        tmpfile.close()
        return disulfides_pairs


    def write_amber_input_file_lib(self,disulfides_pair,libraryfile,parameterfile):
        ds_string = ""
        for i in disulfides_pair:
            pair = i.split(',')
            ds_string = ds_string +"bond prt."+str(pair[0])+".SG prt."+str(pair[1])+".SG\n"
        template = '''
source leaprc.gaff
loadoff '''+libraryfile+'''
loadoff ions08.lib
loadamberparams '''+parameterfile+'''
# load protein
prt = loadpdb design1.pdb
'''+ds_string+'''
saveamberparm prt vacuum.prmtop vacuum.inpcrd
quit'''
        return template

    def write_amber_input_file(self,disulfides_pair):
        ds_string = ""
        for i in disulfides_pair:
            pair = i.split(',')
            ds_string = ds_string +"bond prt."+str(pair[0])+".SG prt."+str(pair[1])+".SG\n"
        template = '''
source leaprc.gaff
loadoff ions08.lib
# load protein
prt = loadpdb design1.pdb
'''+ds_string+'''
saveamberparm prt vacuum.prmtop vacuum.inpcrd
quit'''
        return template



    def write_to_pdb_file(self,pdblist, libraryfile ):
        
        with open('design1.pdb', 'w') as f:
            for line in pdblist:
                # changed 17-06-2014 to type None
                # if( libraryfile == "None" and line[0:4] == "HETA" ):
                if( libraryfile == None and line[0:4] == "HETA" ):
                    continue
                elif( line[0:4] != "ATOM" and line[0:4] != "HETA" ):
                    continue
                else:
                    f.write(line)


    def write_amber_parameterfile(self, libraryfile, parameterfile, sulfide_positions):
        # changed 17-06-2014 to type None
        #        if( libraryfile != "None"):
        if( libraryfile != None ):
            output_file = self.write_amber_input_file_lib( sulfide_positions, libraryfile, parameterfile )
        else:
            output_file = self.write_amber_input_file( sulfide_positions )
            
        with open("generate_input_parameters_disulfides.sh",'w') as f:
            for i in output_file:
                f.write(i)


    def main(self,pdbfile, libraryfile, parameterfile):
        filename = open(pdbfile, 'r')

        tmplist = []
        # HID residues - Ndelta protonation
        hid_residue_id = []
        # HIE residue - Nepsilon protonation
        hie_residue_id = []

        change = False

        for line in filename:
            if ( line[17:20] == "HIS" ):
                # HD1 - proton located on Nd1
                if ( line[13:16] == "HD1" ):
                    hid_residue_id.append(line[23:26].strip())
                # HE2 - proton is located on Ne2
                if ( line[13:16] == "HE2" ):
                    hie_residue_id.append(line[23:26].strip())

            if (line[13:14] != 'H' and line[13:15] != "NV" and line[12:14] != "OV" and line[13:14] != 'X'):
                tmplist.append(line)
        filename.close()
        # replace histidine protonations
        if ( len(hid_residue_id) > 0 or len(hie_residue_id) > 0):
            tmplist = self.set_correct_histidines(tmplist, hie_residue_id)

        # vector with coordinates
        sulfur_coordinates = self.get_sulfur_atom_coordinates(tmplist)

        positions_to_change_for_sulfides = self.get_disulfide_pairs(sulfur_coordinates)

        if ( len(positions_to_change_for_sulfides) > 0):
            tmplist = self.set_disulfide_pairs(tmplist, positions_to_change_for_sulfides)

        self.write_to_pdb_file( tmplist , libraryfile )
        
        self.write_amber_parameterfile(libraryfile, parameterfile, positions_to_change_for_sulfides)

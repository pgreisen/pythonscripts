import os, shutil, sys, subprocess
from numpy import *

'''



'''

DISULFIDESDISTANCE = 2.5

# 
all_atoms = True
names_to_changes = ["XX ", "C2 ", "C3 ", "C4 ", "C5 ", "C6 ", "C7 ", "C8 ", "C9 ", "C10", "O1 ", "H10", "H7 "]

# Input name is key and value is the name from the MD
# parameter file
atom_name = {
    "C3 ": "C  ",
    "O3 ": "O  ",
    "N1 ": "N1 ",
    "N2 ": "N  ",
    "C4 ": "C1 ",
    "C5 ": "C4 ",
    "C2 ": "C2 ",
    "S1 ": "S  ",
    "C6 ": "C3 ",
    "C7 ": "C5 ",
    "C8 ": "C6 ",
    "C9 ": "C7 ",
    "C10": "C8 ",
    "C11": "C9 ",
    "O11": "O2 ",
    "O12": "O1 "
}


def insert_correct_histidines(tmplist, hie_numbers):
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


def insert_disulfides(tmplist, residuepairs):
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



# add 11-02-2014
def get_sulfur_atom_coordinates(pdbfile):
    sulfur_coordinates = {}

    for line in pdbfile:

        if ( line[0:4] == 'ATOM' and line[13:14] == 'S'):
            sulfur_coordinates[line[22:26].strip()] = array(
                [float(line[31:39]), float(line[39:47]), float(line[47:55])])

    return sulfur_coordinates


def test_disulfide(sulfur_coordinates):
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
                if linalg.norm( v - vv ) <= DISULFIDESDISTANCE:
                    tmpfile.write(k+','+kk+"\n")
                    already_done.append(k)
                    pair = k+','+kk

                    #disulfides_pairs.append(k)
                    #disulfides_pairs.append(kk)
                    disulfides_pairs.append(pair)

            else:
                continue

    tmpfile.close()
    return disulfides_pairs


def write_amber_input_file(disulfides_pair):
    # bond oxy.1.SG oxy.6.SG
    ds_string = ""
    for i in disulfides_pair:
        pair = i.split(',')
        ds_string = ds_string +"bond prt."+str(pair[0])+".SG prt."+str(pair[1])+".SG\n"
    # blabla
    template = '''
source leaprc.gaff
#loadoff q4md-forcefieldtools-cofactors.off
#loadoff phosphoaa10.lib
loadoff ions08.lib
# load protein
prt = loadpdb design1.pdb
'''+ds_string+'''

saveamberparm prt apo_vacuum.prmtop apo_vacuum.inpcrd
quit
    '''
    return template

def main():
    pdbfile = sys.argv[1]

    filename = open(pdbfile, 'r')

    tmplist = []
    # HID residues - Ndelta protonation
    hid_residue_id = []
    # HIE residue - Nepsilon protonation
    hie_residue_id = []

    residuename = "BTN"

    new_residue_name = "BTN"

    change = False

    for line in filename:

        if ( line[17:20] == residuename ):

            if (line[13:16] in names_to_changes or all_atoms):

                newline = line[0:13] + str(atom_name[str(line[13:16])]) + line[16:17] + new_residue_name + line[20:]

                tmplist.append(newline)

                # tmpfile.write(newline)

            elif ( line[13:14] == 'H'):

                continue

            elif (line[13:16] in forbidden_atom_names):

                continue

            else:

                newline = line[0:17] + new_residue_name + line[20:]

                tmplist.append(newline)

                # tmpfile.write(newline)

        else:
            # test for histidine

            if ( line[17:20] == "HIS" ):
                # HD1 - proton located on Nd1
                if ( line[13:16] == "HD1" ):
                    hid_residue_id.append(line[23:26].strip())
                    # HE2 - proton is located on Ne2
                if ( line[13:16] == "HE2" ):
                    hie_residue_id.append(line[23:26].strip())

            if (line[13:14] != 'H' and line[13:15] != "NV" and line[12:14] != "OV" and line[13:14] != 'X'):
                tmplist.append(line)


    # tmpfile.close()
    filename.close()

    # replace histidine protonations
    if ( len(hid_residue_id) > 0 or len(hie_residue_id) ):
        tmplist = insert_correct_histidines(tmplist, hie_residue_id)

    # vector with coordinates
    sulfur_coordinates = get_sulfur_atom_coordinates(tmplist)

    positions_to_change_for_sulfides = test_disulfide(sulfur_coordinates)

    if ( len(positions_to_change_for_sulfides) > 0):
        
        tmplist = insert_disulfides(tmplist, positions_to_change_for_sulfides)

    tmpfile = open('tmp.pdb', 'w')

    for line in tmplist:
        tmpfile.write(line)
    tmpfile.close()

    mvfl = 'mv tmp.pdb ' + str(pdbfile)

    subprocess.Popen(mvfl, shell=True).wait()

    tmp_amber = open("generate_input_parameters_disulfides.sh",'w')

    output_file = write_amber_input_file( positions_to_change_for_sulfides )

    for i in output_file:
        tmp_amber.write(i)

    tmp_amber.close()

if __name__ == "__main__":
    main()

import os,shutil,sys,subprocess
'''



'''

def insert_correct_histidines(tmplist, hie_numbers):
    
    correct_list = []

    for line in tmplist:

        if(line[17:19] == "HI" and line[23:26].strip() in hie_numbers):
            
            line = str(line[0:17])+"HIE"+str(line[20:])
        
        elif(line[17:19] == "HI"):

            line = str(line[0:17])+"HID"+str(line[20:])
        
        correct_list.append(line)

    return correct_list



def main():

    pdbfile = sys.argv[1]
    
    filename = open(pdbfile,'r')

    tmplist = []
    # HID residues - Ndelta protonation
    hid_residue_id = []
    # HIE residue - Nepsilon protonation
    hie_residue_id = []
    
    change = False

    for line in filename:


                
        # test for histidine
            
        if( line[17:20] == "HIS" ):
            # HD1 - proton located on Nd1
            if( line[13:16] == "HD1" ):
                hid_residue_id.append( line[23:26].strip() )
            # HE2 - proton is located on Ne2
            if( line[13:16] == "HE2" ):
                hie_residue_id.append( line[23:26].strip() )

        if(line[13:14] != 'H' and line[13:15] != "NV") :
            tmplist.append(line)


    filename.close()

    # replace histidine protonations
    if( len( hid_residue_id ) > 0 or len( hie_residue_id ) ):
        
        correct_list = insert_correct_histidines(tmplist, hie_residue_id)

    tmpfile = open('tmp.pdb','w')
    for line in correct_list:
        tmpfile.write(line)

    tmpfile.close()

    mvfl = 'mv tmp.pdb '+str(pdbfile)

    subprocess.Popen(mvfl,shell=True).wait()

    
if __name__ == "__main__":
    main()

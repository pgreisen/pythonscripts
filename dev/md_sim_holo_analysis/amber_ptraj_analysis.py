__author__ = 'greisen'

import subprocess,os,csv
from numpy import mean,var,sqrt
import math

class amber_ptraj_analysis:

    def __init__(self):
        self.RMSD_CUT_OFF = 1
        self.CHI_ANGLE_SD_CUT_OFF = 20
        self.pdbfile = 0



    # avg_pdb.ptraj
    def average_pdb(self,number_of_trajectories,ligand_number,torsion_angles,rmsd_per_residue,apo):
        print "THE LIGAND NUMBER IS ",ligand_number
        last_protein_residue =  str(int(ligand_number)-1)
        ln = str(ligand_number)
        apo_string = ""
        if(apo != "apo"):
            apo_string = '''
            atomicfluct out rmsf.data :1-'''+str(ln)+''' byres
            atomicfluct out rmsf_ligand.data :'''+str(ln)+''' byres
            diffusion :'''+str(ln)+'''
            '''
        else:
            ln = last_protein_residue

        template = '''
# We use this template
'''+number_of_trajectories+'''
# Remove water
strip :WAT

# first center the solute by itself
center origin :1-'''+str(ln)+'''
# now image the whole system about the centered origin
image origin center

average avg.pdb pdb

quit

        '''


        tmpfile = open("avg_pdb.ptraj",'w')
        for line in template:
            tmpfile.write(line)
        tmpfile.close()


    def analysis_template(self,number_of_trajectories,ligand_number,torsion_angles,rmsd_per_residue,apo):
        self.average_pdb(number_of_trajectories,ligand_number,torsion_angles,rmsd_per_residue,apo)


        print "THE LIGAND NUMBER IS ",ligand_number
        last_protein_residue =  str(int(ligand_number)-1)
        ln = str(ligand_number)
        apo_string = ""
        if(apo != "apo"):
            apo_string = '''
            atomicfluct out rmsf.data :1-'''+str(ln)+''' byres
            atomicfluct out rmsf_ligand.data :'''+str(ln)+''' byres
            diffusion :'''+str(ln)+'''
            # Insert donor and acceptor residues
            donor mask :HCY@O4
            donor mask :HCY@O1

            donor mask :HCY@O2
            donor mask :HCY@O3
            donor mask :HCY@O5

            acceptor mask  :HCY@O5 :HCY@H30
            acceptor mask  :HCY@O3 :HCY@H29
            acceptor mask  :HCY@O2 :HCY@H28

            '''
        else:
            ln = last_protein_residue


        # add the proline string for the HB-analysis to the script
        proline_residues = self.get_proline_positions( self.pdbfile, last_protein_residue )
        print "The proline residues are: ", proline_residues


        template = '''
# We use this template
'''+number_of_trajectories+'''
# Remove water
##strip :WAT

# Removes water, output stripped topology as stripped.prmtop
strip :WAT outprefix stripped

# first center the solute by itself
center origin :1-'''+str(ln)+'''
# now image the whole system about the centered origin
image origin center

##average avg.pdb pdb

reference avg.pdb

'''+torsion_angles+'''
# save the final "new" snapshot as a PDB file
# for alternate veiwing

#trajout last_frame.pdb PDB

# need rms to remove rotation and translation
rms first mass :1-'''+last_protein_residue+'''

# rmsd per residue
'''+rmsd_per_residue+'''
# need to check this command 06-11-2013

'''+apo_string+'''

atomicfluct out rmsf_protein.data :1-'''+str(last_protein_residue)+''' byres

matrix mwcovar out mwcor.dat

#rms reference cc
matrix correl out cross_correlation.dat @CA



# Hydrogen bond analysis
# If a ligand is present you will need
# add those atoms as well for the analysis
#-- Donors from standard amino acids
donor mask :GLN@OE1
donor mask :GLN@NE2
donor mask :ASN@OD1
donor mask :ASN@ND2
donor mask :TYR@OH
donor mask :ASP@OD1
donor mask :ASP@OD2
donor mask :GLU@OE1
donor mask :GLU@OE2
donor mask :SER@OG
donor mask :THR@OG1
donor mask :HIS@ND1
donor mask :HIE@ND1
donor mask :HID@NE2

#-- Acceptors from standard amino acids
acceptor mask  :ASN@ND2 :ASN@HD21
acceptor mask  :ASN@ND2 :ASN@HD22
acceptor mask  :TYR@OH  :TYR@HH
acceptor mask  :GLN@NE2 :GLN@HE21
acceptor mask  :GLN@NE2 :GLN@HE22
acceptor mask  :TRP@NE1 :TRP@HE1
acceptor mask  :LYS@NZ  :LYS@HZ1
acceptor mask  :LYS@NZ  :LYS@HZ2
acceptor mask  :LYS@NZ  :LYS@HZ3
acceptor mask  :SER@OG  :SER@HG
acceptor mask  :THR@OG1 :THR@HG1
acceptor mask  :ARG@NH2 :ARG@HH21
acceptor mask  :ARG@NH2 :ARG@HH22
acceptor mask  :ARG@NH1 :ARG@HH11
acceptor mask  :ARG@NH1 :ARG@HH12
acceptor mask  :ARG@NE  :ARG@HE
acceptor mask  :HIS@NE2 :HIS@HE2
acceptor mask  :HIE@NE2 :HIE@HE2
acceptor mask  :HID@ND1 :HID@HD1
acceptor mask  :HIP@ND1,NE2 :HIP@HE2,HD1

#-- Backbone donors and acceptors for this particular molecule
#   N-H for prolines do not exist so are not in the mask
#
donor mask @O

# Here one has to exclude all the prolines from the analysis.
acceptor mask  :2-11,13-16,20@N :2-20@H
#Terminal residues have different atom names
donor mask @OXT
acceptor mask :1@N :1@H1
acceptor mask :1@N :1@H2
acceptor mask :1@N :1@H3

# -- series hbt is just a placeholder to ensure we get the full analysis. If you don't
# have the word series you don't get a full analysis.
# print only hydrogen bonds that are resident for more than 5% of the simulation
hbond print .05 series hbt




quit

        '''
        return template


    def generate_vmd_input_file(self,residues,ligand_number,number_of_md_simulation=50):
        
        template = '''set start 0
# file name is apo_solv.prmtop
mol load parm7 MD_1/solv.prmtop

set end '''+str(number_of_md_simulation)+'''

for {set start 1} {$start < $end } {incr start} {
    puts "MD_$start"
    mol addfile MD_$start/mdcrd type crdbox waitfor all
}
color Display Background white
# orthographic
display projection orthographic

# Reset the view of the window
mol delrep 0 all
# Make dssp representation of protein
# Numbers at the end are taken from
mol representation NewCartoon
# add color
mol color SecondaryStructure
# The material used for the representation
mol material Opaque
# The selection for all of this
mol selection {all}
# Now update the window
mol addrep top

# show critical residues in CPK
## Color part of selection a certain color
mol representation CPK
# 0 is blue
mol material Opaque
mol selection { residue '''+residues+''' and not hydrogen }
# uncomment to change color
# mol color ColorID 23
mol addrep top

# Load design into trajectory
mol load pdb Minimization/design1.pdb

# align the two
# need to insert this into the analyze script:
# this should align.
# superposition of sel1 onto sel0

set sel0 [atomselect 0 "backbone and resid 1 to '''+str( int(ligand_number)-1 )+''' "]
set sel1 [atomselect 1 "backbone and resid 1 to '''+str( int(ligand_number)-1 )+''' "]
set M [measure fit $sel1 $sel0]

#$sel1 move $M

# move this with an additional residue which is the ligand
# if not the ligand will just stay with its original coordinates.
# set sel2 [atomselect 1 "resid '''+str(ligand_number)+''' "]
set sel2 [atomselect 1 "all "]
$sel2 move $M

mol representation CPK
# 0 is blue
mol material Opaque
mol selection { residue '''+str(ligand_number)+''' }
mol addrep top


        '''
        dump_template_to_file = open("load_files.vmd",'w')
        dump_template_to_file.write( template )
        dump_template_to_file.close()




    # formula used for conversion between B-factor and
    # RMSF RMSF = sqrt(3)*B/(8*pi^2)
    def convert_rmsf_into_bfactor(self):
        bfactor = []
        rmsf = []
        rmsf_file = open("rmsf_protein.data",'r')
        for line in rmsf_file:
            bfactor.append( 8*math.pi**2 *float(line.split()[1])/sqrt(3) )
        return bfactor



    def number_of_trajecties(self,number_of_directories):
        template_string = ""
        assert isinstance(number_of_directories, int)
        i = 0
        while i <= number_of_directories:

            template_string = template_string+"trajin MD_"+str(i)+"/mdcrd\n"
            i += 1
        return template_string

    def get_torsion_analysis(self,protein_position_and_type):
        '''
        dihedral trpchi1 :6@N :6@CA :6@CB :6@CG out trp_chi1.dat
        dihedral trpchi2 :6@CA :6@CB :6@CG :6@CD2 out trp_chi2.dat
        @return:
        '''
        chi_angles_for_amber = ""
        for key in protein_position_and_type:
            chi_angle = 0

            tmp_residue = self.get_chi_angles( protein_position_and_type[key] )
            for torsion_angle in tmp_residue:
                chi_angle += 1
                tmp = ""
                for ta in torsion_angle:
                    tmp = tmp+" :"+str(key)+"@"+str(ta)
                tmp_string = "dihedral "+str(protein_position_and_type[key])+"chi"+str(chi_angle)+"_"+str(key)+" "
                new_string = tmp_string+tmp+" out "+str(key)+"_"+str(protein_position_and_type[key])+"_chi_"+str(chi_angle)+"_chi.dat"
                chi_angles_for_amber = chi_angles_for_amber+ new_string+"\n"

        return chi_angles_for_amber


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
        'HID' : [ ['N','CA','CB','CG'],['CA','CB','CG','ND1'] ],
        'HIE' : [ ['N','CA','CB','CG'],['CA','CB','CG','ND1'] ],
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

    def get_rmsd_per_residue_string(self,protein_position_and_type):
        rmsd_per_residue = ""
        for key in protein_position_and_type:

            tmp_string = "rms reference out "+protein_position_and_type[key]+"_"+str(key)+"_rmsd.dat :"+str(key)+" nofit"

            rmsd_per_residue = rmsd_per_residue+tmp_string+"\n"


        return rmsd_per_residue



    def generate_critical_residues_session(self):

        print "under construction"


    def get_torsion_angles_rosetta_design(self):
        tmp = open("chi_table.txt",'r')
        chi_table = tmp.readlines()
        tmp.close()
        return chi_table


    def get_statistics_rmsd(self,fl):
        tmp = open(fl,'r')
        data = []
        for line in tmp:
            data.append(float(line.split()[1]) )
        filename = fl.split('_')

        return filename[0],filename[1],round(mean(data),2),round(sqrt(var(data)),2),round(min(data),2),round(max(data),2)


    def write_csv_table_rmsd_analysis(self):
        critical_rmsd_residues = []
        # table with chi-values
        with open('rmsd_per_residue_analysis.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["Amino Acid","Residue NR","RMSD","sd","min","max"])
            # rmsd.dat
            path = './'
            files = os.listdir(path)
            for fl in files:
                if os.path.isfile(fl):
                    if( fl.endswith("rmsd.dat") ):
                        resname,residuenr,rmsd,sd,min_v,max_v = self.get_statistics_rmsd(fl)

                        spamwriter.writerow([resname,residuenr,rmsd,sd,min_v,max_v])

                        if( float( max_v ) >  1.0 ): #self.RMSD_CUT_OFF ):

                            tmp_string = "Residue "+resname+" "+residuenr+" with max RMSD "+str(max_v)+"\n"
                            critical_rmsd_residues.append(tmp_string )



        return critical_rmsd_residues


    def get_statistics_torsion(self,fl):
        tmp = open(fl,'r')
        data = []
        for line in tmp:
            # to avoid making the trans-angle
            value = float(line.split()[1])
            if( value <= -140 ):
                value = abs(value)
            if ( value == None):
                break
            data.append(float( value ) )
        filename = fl.split('_')

        return filename[1],filename[0],filename[3],round(mean(data),2),round(sqrt(var(data)),2),round(min(data),2),round(max(data),2)


    def get_difference_torsion(self,resname,residuenr,chi_angle):


        design_chi_table = self.get_torsion_angles_rosetta_design()

        for line in design_chi_table:
            print line, "get_difference_torsion"
            tmp = line.split()

            print tmp[0], resname
            print tmp[1], residuenr

            if( tmp[0] == resname and tmp[1]== residuenr ):
                dummy = int(chi_angle) + 1
                try:
                    design_chi_angle = float( tmp[dummy] )
                    if( design_chi_angle <= -140 ):
                        design_chi_angle = abs(design_chi_angle)
                except:
                    design_chi_angle = "nan"

                return design_chi_angle
            else:
                return "nan"


    def write_csv_table_torsion_analysis(self):
        critical_torsion_residues = []
        # table with chi-values
        with open('chi_analysis.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["Amino Acid","Residue NR","X-anlge","mean","sd","min","max","designed","delta"])
            path = './'
            files = os.listdir(path)
            for fl in files:
                if os.path.isfile(fl):
                    if( fl.endswith("_chi.dat") ):
                        print "chi file is", fl
                        resname,residuenr,chi_angle,mean_v,sd,min_v,max_v = self.get_statistics_torsion(fl)
                        print resname,residuenr,chi_angle,mean_v,sd,min_v,max_v
                        designed_torsion_angle = self.get_difference_torsion(resname,residuenr,chi_angle)
                        print designed_torsion_angle
                        delta = float(designed_torsion_angle) - mean_v
                        if( delta > 20 ): #self.CHI_ANGLE_SD_CUT_OFF):
                            tmp_string = "Residue "+resname+" "+residuenr+" with chi angle "+chi_angle+" with delta "+str(delta)+"\n"
                            critical_torsion_residues.append(tmp_string )
                        spamwriter.writerow([resname,residuenr,chi_angle,mean_v,sd,min_v,max_v,str(designed_torsion_angle),str(delta)])
        return critical_torsion_residues

    def get_number_of_directories(self):
        # loop over all the
        number_of_directories = 0
        path = './'
        dirs = os.listdir(path)
        for dr in dirs:
            if os.path.isdir(dr):
                if( dr.startswith( "MD") ):
                    number_of_directories += 1
        return number_of_directories


    def get_ligand_position(self,pdbfile,ligandname):
        ligandnumber = 0
        for line in pdbfile:
            ligandnumber = str(line[21:26]).strip()

            if( str(line[17:20]).strip() == ligandname ):
                ligandnumber = str( line[21:26].strip())
                return ligandnumber

            elif( str(line[17:20]).strip() == "Na+"): #ligandname ):
                ligandnumber = str( int(line[21:26].strip()) + 1)
                return ligandnumber

            elif( str(line[17:20]).strip() == "Cl-"): #ligandname ):
                ligandnumber = str( int(line[21:26].strip()) + 1)
                return ligandnumber

            elif( str(line[17:20]).strip() == "WAT"): #ligandname ):
                ligandnumber = str( int(line[21:26].strip()) + 1)
                return ligandnumber



        return 0


    def execute_ptraj(self,ligand_parameter):

        generate_avg = "~greisen/amber12/bin/ptraj MD_1/"+str(ligand_parameter)+" < avg_pdb.ptraj"

        # run amber command
        exe = "~greisen/amber12/bin/ptraj MD_1/"+str(ligand_parameter)+" < analysis.ptraj > hbond_analysis.dat"
        subprocess.Popen(generate_avg,shell=True).wait()

        subprocess.Popen(exe,shell=True).wait()


    def write_template_to_file(self,template):
        tmp = open("analysis.ptraj","w")
        for line in template:
            tmp.write(line)


    def get_proline_positions(self,pdbfile, length_of_protein ):
        proline_residues = []
        residue_number = 0

        pair_residues = ''

        for line in pdbfile:

            if line[17:20] == "PRO" and str(line[23:26]).strip() != residue_number :
                residue_number = str(line[23:26]).strip()
                proline_residues.append( residue_number )
        non_proline_residues = []

        i = 1
        while i <= int(length_of_protein):
            if str(i) not in proline_residues:
                non_proline_residues.append(i)
            i = i + 1

        # prune list
        prune_list = []
        pl = True


        for i in range( len( non_proline_residues ) ):

            if(non_proline_residues[i]  == 1 ):
                continue

            elif( pl == True ):
                start_residue = non_proline_residues[i]
                pl = False

            elif( non_proline_residues[i] == non_proline_residues[i-1] + 1 ):
                # print non_proline_residues[i],non_proline_residues[i-1] + 1,non_proline_residues[i] == non_proline_residues[i-1] + 1

                continue

            else:
                pair_residues = ''
                end_residue = non_proline_residues[i-1]
                if( start_residue == end_residue ):
                    pair_residues = start_residue
                else:
                    pair_residues = str(start_residue)+'-'+str(end_residue)
                prune_list.append( pair_residues )
                start_residue = non_proline_residues[i]
        # acceptor mask  :2-11,13-16,20@N :2-20@H
        mask_string = "acceptor mask  :"

        print prune_list, "proline"
        for tmp in prune_list:
            mask_string = mask_string +str(tmp)+','
        if( start_residue != length_of_protein):
            mask_string = mask_string+str(start_residue)+'-'+str(length_of_protein)+'@N :2-'+str(length_of_protein)+"@H"


        return proline_residues

    def read_pdbfile(self):
        tmp = open("tmp.pdb",'r')
        pdbfile = tmp.readlines()
        tmp.close()
        return pdbfile


    def get_parameterfile_and_rst_file(self):

        parameterfile = ""
        rst_file = ""


        # Go into MD_1
        path = './'
        dirs = os.listdir(path)

        for dr in dirs:

            if os.path.isdir(dr):
                if(dr == "MD_1"):

                    os.chdir(dr)

                    files = os.listdir(path)

                    for fl in files:

                        if ( fl.endswith( "prmtop" ) ):
                            parameterfile = fl

                        elif ( fl == "md.rst"):
                            rst_file = fl

                    # execute command and move file
                    exe = "~greisen/amber12/bin/ambpdb -p "+str(parameterfile)+" <"+rst_file+" > tmp.pdb"
                    subprocess.Popen(exe,shell=True).wait()
                    # move file down
                    mv_file = "mv tmp.pdb ../"
                    subprocess.Popen(mv_file,shell=True).wait()
                    os.chdir('../')
                    break
        print "getting pdb file from MD simulation in directory MD_1",parameterfile
        return parameterfile

    def clean_up_directory(self):
        mkdir_data_dir = "mkdir data_files_for_analysis"
        subprocess.Popen(mkdir_data_dir,shell=True).wait()
        mv_data_files = "mv *dat *txt  data_files_for_analysis/"
        subprocess.Popen(mv_data_files,shell=True).wait()
        remove_tmp_pdb = "rm -f tmp.pdb"
        subprocess.Popen(remove_tmp_pdb,shell=True).wait()


    def get_cross_correlation_network(self,residue_number):
        line_number = int(residue_number) - 1
        dummy = 0
        cc_vector = []

        os.chdir("data_files_for_analysis")

        # data_files_for_analysis/cross_correlation.dat
        cc_file = open("cross_correlation.dat",'r')
        for line in cc_file:
            if( dummy == line_number ):
                tmp = line.split()
                for i in tmp:
                    cc_vector.append(i)
            dummy +=1
        os.chdir("../")
        # write to file
        tmp_write = open("cc_vector_"+str(residue_number)+".dat",'w')
        for number in cc_vector:
            tmp_write.write( str(number)+"\n" )
        return cc_vector




    def dump_pymol_file_with_critical_residues(self,torsion_positions,rmsd_positions):
        # insert string with rmsf data
        tp = ""
        cc_tp = []
        for line in torsion_positions:
            tmp_line = line.split()
            if( tmp_line[0] == "Residue"):
                tmp_residue_number = tmp_line[2]
                # getting cc for torsional residues
                cc_tp.append( self.get_cross_correlation_network( tmp_residue_number ) )
                tp = tp+str( tmp_residue_number )+"+"
        cc_rp = []
        rp = ""
        for line in rmsd_positions:
            tmp_line = line.split()
            if( tmp_line[0] == "Residue"):
                tmp_residue_number = tmp_line[2]
                # getting cc for torsional residues
                cc_rp.append( self.get_cross_correlation_network( tmp_residue_number ) )
                rp = rp+str(tmp_residue_number)+"+"

        with open('critical_residues.pml', 'w') as pymol_session:
            pymol_session.write("hide everything\n")
            pymol_session.write("show cartoon\n")
            pymol_session.write("show sticks, het\n")

            pymol_session.write("create protein_rmsf, all and not het\n")

            if( len(tp) > 0):

                pymol_session.write("create torsional_angles, resi "+tp+"\n")
                pymol_session.write("create rmsd_residues, resi "+rp+"\n")

                pymol_session.write("show sticks, torsional_angles\n")
                pymol_session.write("show sticks, rmsd_residues\n")

                pymol_session.write("color red, torsional_angles and name c*\n")
                pymol_session.write("color yellow, rmsd_residues and name c*\n")

            pymol_session.write("inFile = open(\"bfactor_data\", 'r')\n")
            pymol_session.write("stored.newB = []\n")
            pymol_session.write("for line in inFile.readlines(): stored.newB.append( float(line) )\n")
            pymol_session.write("inFile.close()\n")
            pymol_session.write("alter protein_rmsf, b=0.0\n")
            pymol_session.write("alter protein_rmsf and n. CA, b=stored.newB.pop(0)\n")
            pymol_session.write("spectrum b, rainbow_rev, protein_rmsf and n. CA\n")

            pymol_session.write("create p_rmsf, protein_rmsf\n")
            pymol_session.write("spectrum b, rainbow_rev, p_rmsf and n. CA\n")
            pymol_session.write("cartoon putty, p_rmsf\n")

            # cc data
            for cc_rp_residue in rp.split('+'):
                if( cc_rp_residue == ""):
                    continue
                # import pdb; pdb.set_trace()
                pymol_session.write("inFile = open(\"cc_vector_"+str( cc_rp_residue )+".dat\", 'r')\n")
                pymol_session.write("stored.newB = []\n")
                pymol_session.write("for line in inFile.readlines(): stored.newB.append( float(line) )\n")
                pymol_session.write("inFile.close()\n")

                pymol_session.write("create cc_residue_"+str(cc_rp_residue)+", protein_rmsf\n")
                pymol_session.write("alter cc_residue_"+str(cc_rp_residue)+", b=0.0\n")
                pymol_session.write("alter cc_residue_"+str(cc_rp_residue)+" and n. CA, b=stored.newB.pop(0)\n")
                pymol_session.write("spectrum b, rainbow_rev, cc_residue_"+str(cc_rp_residue)+" and n. CA\n")


            pymol_session.write("hide everything, elem h\n")



    def write_critical_positions(self,torsion_positions,rmsd_positions):

        with open('critical_positions.txt', 'w') as critical_file:
            critical_file.write("These rotamers changes: \n")
            for line in torsion_positions:
                critical_file.write(line)
            critical_file.write("RMSD is high for: \n")
            for line in rmsd_positions:
                critical_file.write(line)
        self.dump_pymol_file_with_critical_residues(torsion_positions,rmsd_positions)

        print "Done with whole analysis"

    def dump_bfactor_to_file(self,bfactor):
        tmp = open("bfactor_data",'w')
        for line in bfactor:
            tmp.write(str(line)+"\n")
        tmp.close()



    def get_ptraj_analysis_file(self,ligandname,protein_position_and_type,apo):

        torsion_angles = self.get_torsion_analysis(protein_position_and_type)
        rmsd_per_residue = self.get_rmsd_per_residue_string(protein_position_and_type)

        md_parameter_file = self.get_parameterfile_and_rst_file()
        self.pdbfile = self.read_pdbfile()

        # get ligand number
        ligand_number = self.get_ligand_position(self.pdbfile, ligandname)


        number_of_directories = self.get_number_of_directories()
        string_directories = self.number_of_trajecties(number_of_directories)

        template = self.analysis_template(string_directories,ligand_number,torsion_angles,rmsd_per_residue,apo)

        self.write_template_to_file(template)
        self.execute_ptraj( md_parameter_file )

        critical_torsion_residues = self.write_csv_table_torsion_analysis()
        critical_rmsd_residues = self.write_csv_table_rmsd_analysis()

        # list with bfactor values
        bfactor = self.convert_rmsf_into_bfactor()

        self.dump_bfactor_to_file(bfactor)
        print "Done with analysis"
        self.clean_up_directory()

        self.write_critical_positions(critical_torsion_residues,critical_rmsd_residues)

        # generate a set between the two lists
        critical_positions = list( set( critical_torsion_residues + critical_rmsd_residues ) )
        residues = ''
        for i in critical_rmsd_residues:
            residues = residues + ' '+str(int(i.split()[2]) -1)

        self.generate_vmd_input_file(residues,ligand_number)



import os,shutil,sys
from align_to_sulfor import *

class set_alignment:
# Definition of treshold for rmsd it is set to 0.5

	
	"""""
	All input paramters are defined here.
	
	"""""
	
	# Requires floating point number
	# Returns floating point with correct
	# number of digits for pdb
	def set_number_digits(self,number):	
		return '%.3f' %number

	def set_length_digit(self,number):
		lngth = len(number)
		if lngth == 7:
			return ' '+number
		if lngth == 6:
			return '  '+number
		if lngth == 5:
			return '   '+number
		if lngth == 4:
			return '    '+number
		else:
			return number

	# Setting path to substrate
	# Requires value 0 - 2 which corresonds to the different substrates at the moment
	# Returns path to files and correct pdbfile
	def set_path_substrate(self,value):
		# 0
		if value == 0:
			# Contains the path to the rotamer enzemble
			path = '/Users/greisen/Data/EnzymeDesign/NewProteins/BinuclearPDBs/Amiton/TS_structures/Rotamers/'
			# The parametrised pdb file
			rosetta = '/Users/greisen/Data/EnzymeDesign/NewProteins/BinuclearPDBs/Amiton/TS_structures/LG_0001.pdb'
		if value == 1:
			path = '/Users/greisen/Substrates/Binuclear/Hydroxy_4atoms/Rotamers/'
			rosetta = '/Users/greisen/Substrates/Binuclear/Hydroxy_4atoms/LG_0001.pdb'
		if value == 2:
			path = '/Users/greisen/Substrates/Binuclear/PhosphorylAlignment/Rotamers/'
			rosetta = '/Users/greisen/Substrates/Binuclear/PhosphorylAlignment/LG_0001.pdb'
		if value == 3:
			path = '/Users/greisen/Substrates/Binuclear/Phosphoryl4Alignment/Rotamers/'
			rosetta = '/Users/greisen/Substrates/Binuclear/Phosphoryl4Alignment/LG_0001.pdb'		
		return path,rosetta

	# Requires list of ligand atoms  
	# Writes file 'tmp.pdb' for the alignment
	def set_align_file(self,atm,ALIGNMENT):
		# Same geometry is used for all alignment
		if(ALIGNMENT == 0):
			pth = '/Users/greisen/Substrates/Binuclear/HydroxyAlignment/hydroxy_alignment.pdb'
		if(ALIGNMENT == 1):
			pth = '/Users/greisen/Substrates/Binuclear/Hydroxy_4atoms/hydroxy_alignment.pdb'
		if(ALIGNMENT == 2):
			pth = '/Users/greisen/Substrates/Binuclear/HydroxyAlignment/hydroxy_alignment.pdb'
		if(ALIGNMENT == 3):
			pth = '/Users/greisen/Substrates/Binuclear/Phosphoryl4Alignment/phosphoryl4alignment.pdb'
		fileobject = open(pth,'r')
		# Get atoms to use for the alignment
		# Written to temporary file
		fname = 'tmp.pdb'
		tm = open(fname,'w')
		for line in fileobject:
			if len(line.split()) > 3:
				atm_name = line.split()[2]
				if atm_name in atm:
					tm.write(line)
		tm.close()

	# took directory with conformations of ligand ensemble
	# if we generate the ensemble after the alignment this is not
	# necessary
	# Returns a list with transformed coordinates
	def get_aligned_coor(self,path,VIZ,templateFile,crystal_coor):
		RMSD_TRESHOLD = 0.3
		obj = align_to_substrate()
		files = os.listdir(path)
		# Reading data from crystal structure where one wants
		# the alignment from
		cry_data,atom_names = obj.get_data(crystal_coor)
		###############################################
		# som vi har lavet for vores model system
		fname = 'tmp.pdb'
		# her faar vi navne 
		sub_data,atom_names = obj.get_data(fname)
		# Superimpose substrate data in crystal structure
		# getting the translation and rotation matrix
		t_m, r_m = obj.get_rotate_translate(sub_data,cry_data)
		# Getting the transformed coordinates
		nw = obj.get_transformed_coor(sub_data,cry_data)
		# We transform the original data
		rmsd_align = obj.get_rmsd(nw,cry_data)
		print 'RMSD score of alignment',rmsd_align
		# Name of list the transformed coordinates are collected in
		outfile = []
		if rmsd_align < RMSD_TRESHOLD:			
			for fl in files:
				ph = path+'/'+fl
				rd = open(ph,'r')
				# Hvad indeholder denne file og hvor er den genereret
				# Filen indeholder data fra modellen, altsaa de data
				# som vi har lavet for vores model system
				fname = 'tmp.pdb'
				# her faar vi navne 
				sub_data,atom_names = obj.get_data(fname)
				# Superimpose substrate data in crystal structure
				# getting the translation and rotation matrix
				t_m, r_m = obj.get_rotate_translate(sub_data,cry_data)
				# Getting the transformed coordinates
				nw = obj.get_transformed_coor(sub_data,cry_data)
				# We transform the original data
				sub,at = obj.get_data(ph)
				# The transformed coordinates
				# what is their construction
				t_c = dot(sub,r_m)+t_m
				# Writing the coordinates
				# Files name of coordinates is
				# Writing to a file called superimposed.pdb
				obj.write_pdb(at,t_c)
				# Rosetta naming convention file which is generated
				# earlier.
				# File for rosetta with the correct naming
				# I/O of file
				sp_file = open('superimposed.pdb','r')
				rosetta = open(templateFile,'r')
				fileOne = sp_file.readlines()
				fileTwo = rosetta.readlines()
				rosetta.close()
				# Variable to count line number in other file
				# used to insert at the right line
				ct = 0
				for i in fileTwo:
					ln = fileOne[ct].split()
					# A very temporary fix for the number greater 100
					x = self.set_number_digits(float(ln[6]))
					y = self.set_number_digits(float(ln[7]))
					z = self.set_number_digits(float(ln[8]))
					x = self.set_length_digit(x)
					y = self.set_length_digit(y)
					z = self.set_length_digit(z)		
					i = str(i[0:30])+x+y+z+str(i[55:81])
					outfile.append(i)
					ct = ct +1
				outfile.append(VIZ)
		else:
			outfile = 'RMSD TOO HIGH'
		return outfile


def main(VZ=1):
    # Parameter to set substrate structures
    # 0: p-nitrophenyl isomer 1
    # 1: p-nitrophenyl isomer 2
    # 2: carbaryl 
    substrate =  int(sys.argv[1])
    
    # Atoms used for the alignment
    # Give them as string names
    atm = [] 
    for arg in sys.argv[2:]:
        atm.append(str(arg))
    # String to determine the format
    # 1. for rosetta  2. VMD
    if VZ == 1:
        VIZ = 'ENDMLD\n model\n'
    else:
        VIZ = 'TER\n'

    # Path to substrate coordinates
    path,templateFile = set_path_substrate(substrate)
    set_align_file(substrate,atm)
    get_aligned_coor(path,atm,VIZ,templateFile)



if __name__ == "__main__":
    main()


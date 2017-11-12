import os,shutil,sys
from align_to_sulfor import *

class set_alignment:

	
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


	# How to optimize this part of the code???
	# Method to get data from each rotamer

	def get_data_to_align(self,filename,SUBSTRATE_ATOMS):
		dic = {}
		for line in filename:
			tmp = line.split()
			if tmp[2] in SUBSTRATE_ATOMS:
				dic[tmp[2]] = line

		wrt = open('tmp.pdb','w')
		
		for i in SUBSTRATE_ATOMS:
			wrt.write(str(dic[i]))

		wrt.close()


	# took directory with conformations of ligand ensemble
	# if we generate the ensemble after the alignment this is not
	# necessary
	# Returns a list with transformed coordinates
	def get_aligned_coor(self,VIZ,LIGANDFILE,crystal_coor,RMSD_TRESHOLD,SUBSTRATE_ATOMS):
		RMSD_BOOLEAN = 0

		obj = align_to_substrate()
		#		files = os.listdir(path)
		outfile = []
		# Reading data from crystal structure where one wants
		# the alignment from
		cry_data,atom_names = obj.get_data(crystal_coor)


		rd = open(LIGANDFILE,'r')
		
		# Hvad indeholder denne file og hvor er den genereret
		# Filen indeholder data fra modellen, altsaa de data
		# som vi har lavet for vores model system

		self.get_data_to_align(rd,SUBSTRATE_ATOMS)
		fname = 'tmp.pdb'
		# her faar vi navne 
		sub_data,atom_names = obj.get_data(fname)
		# Superimpose substrate data in crystal structure
		# getting the translation and rotation matrix
		t_m, r_m = obj.get_rotate_translate(sub_data,cry_data)
		# Getting the transformed coordinates
		nw = obj.get_transformed_coor(sub_data,cry_data)
		rmsd_align = obj.get_rmsd(nw,cry_data)

		
		if rmsd_align < RMSD_TRESHOLD:
			RMSD_BOOLEAN = 1
			sub,at = obj.get_data(LIGANDFILE)

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
			rosetta = open(LIGANDFILE,'r')
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
		return RMSD_BOOLEAN,outfile

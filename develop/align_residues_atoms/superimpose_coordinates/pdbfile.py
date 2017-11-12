import os,shutil,sys
from translate_rotate import *

class pdbfile:

	'''

	read_file : readpdb file and returns list
	get_atoms    : takes a pdb list and specific atoms and returns them 
	replace_xyz  : takes a file and transform its coordinates
	
	'''
	

	def write_file(self,list_of_lines,filename='new.pdb'):
		fl = open(filename,'w')
		for line in list_of_lines:
			fl.write(line)
		fl.close()
		




	# Requires floating point number
	# Returns floating point with correct
	# number of digits for pdb
	def set_number_digits(self,number):	
		return '%.3f' %number


	# @Requires pdbfile
	# @Returns file as list
	def read_file(self,file):
		f = open(file,'r')
		fl = f.readlines()
		f.close()
		return fl

	# @Requires pdblist
	# @Requires list with the desired atoms 
	# @Requires residue - necessary to specify residue name and number
	def get_atoms(self,pdblist,list_of_atomnames,residuename):

		dict_of_coordinates = {}
		
		for line in pdblist:

			if str(line[17:20]).strip() == residuename:
				atomname = str(line[12:16]).strip() 
				
				if atomname in list_of_atomnames:
				
					coordinates = self.get_xyz(line)

					dict_of_coordinates[atomname] = coordinates

		return dict_of_coordinates


	# @Requires pdblist
	# @Requires list with the desired atoms 
	# @Requires residue - necessary to specify residue name and number
	def get_all_residueatoms(self,pdblist,residuename):
	
		dict_of_coordinates = {}
		
		for line in pdblist:

			if str(line[17:20]).strip() == residuename:

				atomname = str(line[12:16]).strip() 
				coordinates = self.get_xyz(line)
				dict_of_coordinates[atomname] = coordinates

		return dict_of_coordinates


	def get_protein_coordinates(self,pdblist):
		protein = []
		for line in pdblist:
			if line[0:4] == 'ATOM':
				protein.append(line)
		return protein


	# Requires a number
	# Set the right length for the number
	def set_length_digit(self,number):
		number = str(number)
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

	
	# Require file
	# Input crystal coordinates
	def get_transformed_coordinates(self,coordinates,pdbfile):
		t_coordinates = []
		for line in pdbfile:

			atomname = str(line[12:16]).strip()

			if atomname in coordinates:

				x = '%.2f' %coordinates[atomname][0]
				y = '%.2f' %coordinates[atomname][1]
				z = '%.2f' %coordinates[atomname][2]

				x =  self.set_length_digit(x)
				y =  self.set_length_digit(y)
				z =  self.set_length_digit(z)

				line = str(line[0:30])+x+y+z+str(line[55:])
			t_coordinates.append(line)
		return t_coordinates


        def get_xyz(self,pdbline):
                x = float(pdbline[30:39])
                y = float(pdbline[39:47])
                z = float(pdbline[47:55])

		return array([x,y,z])


		

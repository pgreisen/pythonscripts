import os,shutil,sys
from translate_rotate import *

class pdbfile:

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

	# Method to get data from each rotamer
	def get_data_to_align(self,filename):
		tmp_chi = open(filename,'r')

		atoms = ['ZN1','ZN2','O5','O1']
		dic = {}
		for line in tmp_chi:
			tmp = line.split()
			if tmp[2] in atoms:
				dic[tmp[2]] = line
		wrt = open('tmp.pdb','w')

		wrt.write(str(dic['ZN2']))
		wrt.write(str(dic['ZN1']))
		wrt.write(str(dic['O5']))
		wrt.write(str(dic['O1']))
		wrt.close()


	# took directory with conformations of ligand ensemble
	# if we generate the ensemble after the alignment this is not
	# necessary
	# Returns a list with transformed coordinates
	def get_aligned_coor(self,path,VIZ,templateFile,crystal_coor):
		RMSD_TRESHOLD = 0.8
		obj = align_to_substrate()
		files = os.listdir(path)
		outfile = []
		# Reading data from crystal structure where one wants
		# the alignment from
		cry_data,atom_names = obj.get_data(crystal_coor)
		for fl in files:
			ph = path+'/'+fl
			rd = open(ph,'r')
			# Hvad indeholder denne file og hvor er den genereret
			# Filen indeholder data fra modellen, altsaa de data
			# som vi har lavet for vores model system
			self.get_data_to_align(ph)
			fname = 'tmp.pdb'
			# her faar vi navne 
			sub_data,atom_names = obj.get_data(fname)
			# Superimpose substrate data in crystal structure
			# getting the translation and rotation matrix
			t_m, r_m = obj.get_rotate_translate(sub_data,cry_data)
			# Getting the transformed coordinates
			nw = obj.get_transformed_coor(sub_data,cry_data)
			rmsd_align = obj.get_rmsd(nw,cry_data)
			print 'rmsd_align',rmsd_align
			print 'rmsd ', rmsd_align
			if rmsd_align < RMSD_TRESHOLD:
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
		return outfile

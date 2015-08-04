# -*- coding: utf-8 -*-
import pymol
from pymol import cmd
import random
 
def findSurfaceResidues(objSel="(all)", cutoff=2.5, doShow=False, verbose=False):
	"""
	findSurfaceResidues
		finds those residues on the surface of a protein
		that have at least 'cutoff' exposed A**2 surface area.
 
	PARAMS
		objSel (string)
			the object or selection in which to find
			exposed residues
			DEFAULT: (all)
 
		cutoff (float)
			your cutoff of what is exposed or not. 
			DEFAULT: 2.5 Ang**2
 
		asSel (boolean)
			make a selection out of the residues found
 
	RETURNS
		(list: (chain, resv ) )
			A Python list of residue numbers corresponding
			to those residues w/more exposure than the cutoff.
 
	"""
	tmpObj="__tmp"
	cmd.create( tmpObj, objSel + " and polymer");
	if verbose!=False:
		print "WARNING: I'm setting dot_solvent.  You may not care for this."
	cmd.set("dot_solvent");
	cmd.get_area(selection=tmpObj, load_b=1)
 
	# threshold on what one considers an "exposed" atom (in A**2):
	cmd.remove( tmpObj + " and b < " + str(cutoff) )
 
	stored.tmp_dict = {}
	cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
	exposed = stored.tmp_dict.keys()
	exposed.sort()
 
        randstr = str(random.randint(0,10000))
	selName = "exposed_atm_" + randstr
	if verbose!=False:
		print "Exposed residues are selected in: " + selName
	cmd.select(selName, objSel + " in " + tmpObj ) 
        selNameRes = "exposed_res_" + randstr
        cmd.select(selNameRes, "byres " + selName )
 
 
	if doShow!=False:
		cmd.show_as("spheres", objSel + " and poly")
		cmd.color("white", objSel)
		cmd.color("red", selName)
 
	cmd.delete(tmpObj)
 
	return exposed
 
 
cmd.extend("findSurfaceResidues", findSurfaceResidues)

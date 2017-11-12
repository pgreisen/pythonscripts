import sys,os


import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
from pymol import cmd,stored

pymol.finish_launching()

# positions verified
exposed = []

'''
tmpObj="__tmp"
cmd.create( tmpObj, objSel + " and polymer");
if verbose!=False:
print "WARNING: I'm setting dot_solvent.  You may not care for this."
cmd.set("dot_solvent");
cmd.get_area(selection=tmpObj, load_b=1)
 
# threshold on what one considers an "exposed" atom (in A**2):
cmd.remove( tmpObj + " and b < " + str(cutoff) )

'''


# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]


'''
# use solvent-accessible surface with high sampling density
set dot_solvent, 1
set dot_density, 3
 
# measure the components individually storing the results for later
alpha1_area=cmd.get_area("alpha1")
beta1_area=cmd.get_area("beta1")
 
 
# measure the alpha1,beta1 pair
ab1_area=cmd.get_area("ab1")
 
# now print results and do some maths to get the buried surface
print alpha1_area
print beta1_area
print ab1_area
print (alpha1_area + beta1_area) - ab1_area
'''

pymol.cmd.load(spath, sname)

pymol.cmd.set("dot_solvent",1);
pymol.cmd.set("dot_density",3);

protein="all and not het"
ligand="het"
complex="all"
ligand_area=pymol.cmd.get_area(selection=ligand, load_b=1)
protein_area=pymol.cmd.get_area(selection=protein, load_b=1)
complex_area=pymol.cmd.get_area(selection=complex, load_b=1)

##print ligand_area
##print protein_area
##print complex_area

print (ligand_area + protein_area ) - complex_area

pymol.cmd.quit()

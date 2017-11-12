import sys,os


import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
from pymol import cmd,stored

pymol.finish_launching()

# changed how exposed the residue should
lower_threshold= 25

# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]
pymol.cmd.load(spath, sname)
pymol.cmd.set("dot_solvent");
tmpObj="all and not elem h and not name c+o+n+ca"
pymol.cmd.get_area(selection=tmpObj, load_b=1)
sasa_residues = cmd.remove( tmpObj + " and b < " + str(lower_threshold))

stored.m = [] ;
cmd.iterate(tmpObj,"print b,name,resn,resi" )
cmd.iterate(sasa_residues,"print b,name,resn,resi" )
cmd.iterate(tmpObj,'stored.m.append(resi)' )

with open("exposed_residues.pos",'w') as f:
    for residue in set( stored.m ):
        f.write( str(residue)+" ")

pymol.cmd.quit()

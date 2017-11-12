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

def get_positions(positionfile):
    positions = ""
    with open(positionfile,'r') as f:
        for line in f:
            
            tmp_line = line.split()
            for character in tmp_line:
                positions = positions + character+"+"
    return positions


# threshold
threshold=4.5

# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]

positionfile = os.path.abspath(sys.argv[2])

positions = get_positions(positionfile)


pymol.cmd.load(spath, sname)

pymol.cmd.set("dot_solvent");

tmpObj="resi "+positions

pymol.cmd.get_area(selection=tmpObj, load_b=1)

# threshold on what one considers an "exposed" atom (in A**2):
sasa_residues = cmd.remove( tmpObj + " and b < " + str(threshold) )

# print "Exposed ",len(sasa_residues)

stored.tmp_dict = {}
cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
exposed = stored.tmp_dict.keys()
exposed.sort()

print exposed
# write to file
with open("exposed_residues.pos",'w') as f:
    for pair in exposed:
        print pair[1]
        f.write( str(pair[1])+" ")


##pymol.cmd.disable("all")
##pymol.cmd.enable(sname)
##pymol.cmd.png("my_image.png")
# Get out!
pymol.cmd.quit()

import sys,os
import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
from pymol import cmd,stored

pymol.finish_launching()
# positions verified
exposed = []
# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]
pymol.cmd.load(spath, sname)
tmpObj="all and name ca and ss s"
stored.tmp_dict = {}
cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
exposed = stored.tmp_dict.keys()
exposed.sort()
for key in exposed:
    tmp_dis= cmd.distance("name CA and resi "+str(key[1])," resn TSV and name O8 ")
    if(tmp_dis < 10):
        try:
            print spath.split('/')[-1], tmp_dis
            break
        except:
            continue
# Get out!
pymol.cmd.quit()

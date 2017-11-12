#!/usr/bin/env python
import sys, shutil, os, subprocess

pdbname = sys.argv[1]

pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+pdbname+'_nohet_1.pdb'

move_files = 'scp dig1:'+pdb_path+' .'

subprocess.Popen(move_files,shell=True).wait()



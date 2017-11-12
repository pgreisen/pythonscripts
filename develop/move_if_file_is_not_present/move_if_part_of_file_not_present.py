#!/usr/bin/env python
import sys, shutil, os, subprocess

'''
Move a directory if a file is present

'''

path = './'

missing_files = open('missingjobs.txt','w')

not_submitted = open('not_submitted.txt','w')

files = os.listdir(path)

file_present = sys.argv[1]


try:
    os.mkdir('Initial_pdbs')
except:
    print "Directory already exists"

for fl in files:

    check = 1

    if os.path.isdir(fl):

        os.chdir(fl)
        
        files2 = os.listdir(path)

        for fl2 in files2:

            if fl2.endswith(file_present):

                check = 0

        os.chdir('../')

        # Not file for running located
        if(check == 0):
            not_submitted.write(fl+"\n")
            # subprocess.Popen('condor_submit '+fl+'/condor.submit',shell=True)

    else:
        continue

#!/usr/bin/env python
import sys, shutil, os, subprocess

NUMBER_OF_FILES = 10000

i = 0

files = os.listdir("./")

for file in files:
    if(i%NUMBER_OF_FILES == 0):
        dst = "dir_"+str(i)
        os.mkdir(dst)
    shutil.move(file,dst+"/")
    i = i + 1

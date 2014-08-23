#!/bin/bash
path=~/rosetta/
nohup $path/rosetta_source/bin/EnzdesFixBB.linuxiccrelease -ex1 -ex2 @enzdes.flags -in:file:s $1 -database $path/modified_rosetta_database -relax:cartesian -dun10 > log &

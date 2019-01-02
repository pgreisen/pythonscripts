#!/bin/bash
#PBS -N BASIC
#PBS -l nodes=1:ppn=1
#PBS -l vmem=1gb
#PBS -l walltime=1:00:00
#PBS -j oe
cd $PBS_O_WORKDIR
~/bin/rosetta_bin/rosetta_scripts.static.linuxgccrelease -database ~/database @$path/flags -parser:protocol $path/rosetta_cm.xml -out:prefix $prefix

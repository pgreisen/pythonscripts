#!/bin/bash
#PBS -N BASIC
#PBS -l nodes=1:ppn=1
#PBS -l vmem=1gb
#PBS -l walltime=1:00:00
#PBS -j oe
# resid_chains format e.g. 10A or 42B
# resname e.g SER
rosetta_pth=~/;
rosetta_db=~/database;
$rosetta_pth/rosetta_scripts.static.linuxgccrelease \
    -database $rosetta_db -parser:protocol $xml @flags \
    -in:file:s $native -out:prefix $prefix \
    -nstruct 1 -out:file:silent_struct_type binary -out:file:silent silent.out \
    -out:file:scorefile -out:prefix $prefix

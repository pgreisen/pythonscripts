#!/bin/bash
#PBS -N BASIC
#PBS -l nodes=1:ppn=1
#PBS -l vmem=1gb
#PBS -l walltime=1:00:00
#PBS -j oe
rosetta_pth=~/;
rosetta_db=~/database;
$rosetta_pth/rosetta_scripts.static.linuxgccrelease -database $rosetta_db -extra_res_fa $parameters -parser:protocol dock.xml -in:file:s $native -in:file:native $native -out:prefix $prefix -nstruct 50 -out:file:silent_struct_type binary -out:file:silent $prefix\_silent.out -out:file:scorefile $prefix\_dock.sc -packing:ex1 -packing:ex2 -out:prefix $prefix @flags -parser:script_vars x=$x -parser:script_vars y=$y -parser:script_vars z=$z
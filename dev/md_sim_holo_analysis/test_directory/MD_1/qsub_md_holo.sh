#!/bin/bash
## RENAME FOR YOUR JOB
#PBS -N MinimizationSolvate
## EDIT FOR YOUR JOB
#PBS -l walltime=4:00:00
## For 12 core nodes are used by the Baker lab
#PBS -l nodes=1:ppn=12,mem=22gb,feature=12core
## EDIT FOR YOUR JOB
## Put the STDOUT and STDERR from jobs into the below directory
## 
## Put both the stderr and stdout into a single file
#PBS -j oe
#PBS -M horst.lechner@uni-graz.at
## a mail is sent when the job is aborted by the batch system.
## b mail is sent when the job begins execution.
## e mail is sent when the job terminates.
## PBS -m abe

## Output directory

## EDIT FOR YOUR JOB
## Sepcify the working directory for this job bundle
## NOTE: we run from the LINKS directory - NOT the dir
## containing the actual script files.  This is CRITICAL
## for our checkpoinitng scheme
#PBS -d ./

mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i md.in -p holo_solv.prmtop -c heat.rst -r md.rst
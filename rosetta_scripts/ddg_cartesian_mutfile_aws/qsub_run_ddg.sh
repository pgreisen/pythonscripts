#!/bin/bash -f
#PBS -N JobName
#PBS -l walltime=72:00:00
#PBS -l mem=10gb
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
###########################################################
cd $PBS_O_WORKDIR;
~/cartesian_ddg.static.linuxgccrelease -database ~/database -s ../lowest.pdb -ddg:mut_file $PBS_O_WORKDIR/mutfile -ddg:iterations 10 -ddg::cartesian -ddg::dump_pdbs false -ddg:bbnbr 1 -fa_max_dis 9.0 -score:weights ref2015_cart -out:path:score $PBS_O_WORKDIR -unmute all
#!/bin/bash -f
#PBS -N JobName
#PBS -l walltime=72:00:00
#PBS -l mem=10gb
~/cartesian_ddg.static.linuxgccrelease -database ~/database -s ./lowest.pdb -ddg:mut_file ./mutfile -ddg:iterations 10 -ddg::cartesian -ddg::dump_pdbs false -ddg:bbnbr 1 -fa_max_dis 9.0 -score:weights ref2015_cart -out:path:score $pth -unmute all
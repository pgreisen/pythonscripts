#!/bin/bash -f
~/cartesian_ddg.static.linuxgccrelease -database ~/database -s $pth/lowest.pdb -ddg:mut_file $pth/$mutfile -ddg:iterations 20 -ddg::cartesian -ddg::dump_pdbs false -fa_max_dis 9.0 -score:weights ref2015_cart -out:path:score $pth -unmute all -ex1 -ex2

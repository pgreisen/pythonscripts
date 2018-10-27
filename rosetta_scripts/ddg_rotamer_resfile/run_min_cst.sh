#!/bin/bash
exe=~/minimize_with_cst.static.linuxgccrelease;
database=~/database/;
$exe -in:file:l lst  -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -database $database -ddg::harmonic_ca_tether 0.5 -score:weights ref2015 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false > mincst.log;
echo "Done with cst minimization";
grep ^c-alpha mincst.log | awk '{print "AtomPair CA "$6" CA "$8" HARMONIC "$10" "$13}' > input.cst;
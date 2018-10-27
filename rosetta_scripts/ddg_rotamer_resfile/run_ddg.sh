exe=~/ddg_monomer.static.linuxgccrelease;
database=~/database/;
inputpdb=../min_w_cst/min_cst_0.5.1_0001.pdb;
inputcst=../min_w_cst/input.cst;
$exe -database $database -resfile resfile -in:file:s $inputpdb -ddg:minimization_scorefunction ref2015 -ddg::iterations 50 -ddg::dump_pdbs false -ignore_unrecognized_res -ddg::local_opt_only false -ddg::min_cst true -constraints::cst_file $inputcst -in::file::fullatom -ddg::min true -ddg::sc_min_only false -ddg::ramp_repulsive true -unmute core.optimization.LineMinimizer -ddg::output_silent true
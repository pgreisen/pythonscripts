#!/bin/bash
~/relax.static.linuxgccrelease -s 1.pdb -use_input_sc \
-constrain_relax_to_start_coords -ignore_unrecognized_res \
-nstruct 20 \
-relax:coord_constrain_sidechains  \
-relax:cartesian -score:weights ref2015_cart \
-relax:min_type lbfgs_armijo_nonmonotone \
-relax:script cart2.script -database ~/database
#!/bin/bash
~/Programs/Rosetta/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -database ~/Programs/Rosetta/Rosetta/main/database/ -parser:protocol mutate_residues.xml @flags -in:file:s $1 -nstruct 1 

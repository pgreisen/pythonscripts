path=/home/ubuntu/LigandPSSM_v3/designfiles
rosetta_pth=/home/ubuntu;
rosetta_db=/home/ubuntu/database;
	$rosetta_pth/rosetta_scripts.static.linuxgccrelease -database $rosetta_db -extra_res_fa $1 -parser:protocol $path/liganddesign_pssm.xml -in:file:s $2 -in:file:native $2 -nstruct 20 -out:file:scorefile ${2%.pdb}.sc -packing:ex1 -packing:ex2 -out:prefix $3 @$path/flags

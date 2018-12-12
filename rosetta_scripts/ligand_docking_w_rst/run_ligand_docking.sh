rosetta_pth=~/;
rosetta_db=~/database;
$rosetta_pth/rosetta_scripts.static.linuxgccrelease -database $rosetta_db -extra_res_fa $1 -parser:protocol dock.xml -in:file:s $2 -in:file:native $2 -out:prefix $3 -nstruct 25 -out:file:silent_struct_type binary -out:file:silent $3\_silent.out -out:file:scorefile $3\_dock.sc -packing:ex1 -packing:ex2 -out:prefix $3 @flags

path_to_pdb=/z/insulin/users/pjug/databases/pdb/*;

# loop over all the subdirectories in the database

for j in $path_to_pdb;
do
    cd $j;
    for i in *gz; #pdb;
    do
	echo $i;
	# ~/Programs/tmalign/TMalign /z/insulin/users/pjug/Projects/factors/TrypsinFamily/20151207_HumanStructures/ProteaseDomainchainA/ChainH/1DAN.pdbH.pdb $i > $i\_logfile;
    done
    cd ..;
done
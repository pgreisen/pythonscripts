


for j in Chain{A..Z};
do
    cd $j;
    for i in *pdb;
    do
	~/Programs/tmalign/TMalign /z/insulin/users/pjug/Projects/factors/TrypsinFamily/20151207_HumanStructures/ProteaseDomainchainA/ChainH/1DAN.pdbH.pdb $i > $i\_logfile;
    done
    cd ..;
done
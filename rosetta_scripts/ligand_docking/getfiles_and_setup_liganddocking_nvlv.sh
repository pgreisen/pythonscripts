number_of_trajectories=25;
name_exe=rosetta_scripts.static.linuxgccrelease;
# Get rosetta executables and database
pth_on_S3=s3://enevolvcomputationalbiology/programs;
pth_on_S3_pdbs=s3://tempfilesdocking/;
pth=`pwd`
if [ ! -d "database" ]; then
    aws s3 cp $pth_on_S3/database.tgz .;
    tar zxf database.tgz;
    database=$pth/database;
    rm database.tgz;
fi

if [ ! -f $name_exe ]; then
    aws s3 cp $pth_on_S3/$name_exe.tgz . ;
    tar zxf $name_exe.tgz;
    rm $name_exe.tgz;
    exe=$pth/$name_exe;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/dock.xml .;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/flags .;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/qsub.sh .;
fi


# assume file contains all necessary files for 
# docking
aws s3 cp $pth_on_S3_pdbs . --include="*.zip" --recursive;
# generalize it
for tmp in *zip;
do
    unzip $tmp;
    rm $tmp;
done

parameters=`find *fa.params`;
x=`sed -n '1p' start_from.txt`;
y=`sed -n '2p' start_from.txt`;
z=`sed -n '3p' start_from.txt`;
pth=`pwd`;

for pdb in *.pdb;
do
    native=$pdb;
    for i in $(seq 1 $number_of_trajectories);
    do
	qsub -v parameters="$pth/$parameters",prefix="$i",x="$x",y="$y",z="$z",native="$pth/$native" qsub.sh;
    done
    echo "Done with $native";
done
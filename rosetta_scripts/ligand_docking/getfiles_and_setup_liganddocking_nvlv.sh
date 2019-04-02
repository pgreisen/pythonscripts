###############################################
hhsearchdir=$1;
aws s3 cp s3://tempfilesdocking/$hhsearchdir.zip . 
unzip $hhsearchdir.zip;
rm $hhsearchdir.zip;
###############################################
aws s3 cp s3://enevolvcomputationalbiology/programs/rosetta_scripts.static.linuxgccrelease.tgz .
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
aws s3 cp s3://enevolvcomputationalbiology/programs/database.tgz .;
tar zxf database.tgz;
rm database.tgz;
#########################################
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/qsub.sh;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/run_setup.sh;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/flags;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_docking/dock.xml;

cd $hhsearchdir;
parameters=`find *params`;
native=`find *pdb`;
x=`sed -n '1p' xyz.txt`;
y=`sed -n '2p' xyz.txt`;
z=`sed -n '3p' xyz.txt`;



for i in $(seq 1 100);
do
    qsub -v parameters=$parameters,prefix=$i,x=$x,y=$y,z=$z,native=$native ../qsub.sh;
done
cd ..;
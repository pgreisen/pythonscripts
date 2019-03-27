###############################################
sudo apt-get update -y
sudo apt install emacs25 -y
###############################################
# this should be the name of the pdb file that you want to 
# compute on
disulfide=$1;
aws s3 cp s3://temppdb/$disulfide.zip . 
unzip $disulfide.zip;
rm $disulfide.zip;
###############################################
mkdir rosetta_bin;
cd rosetta_bin;
aws s3 cp s3://enevolvcomputationalbiology/programs/rosetta_scripts.static.linuxgccrelease.tgz .
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
cd ../..;
aws s3 cp s3://enevolvcomputationalbiology/programs/database.tgz .;
tar zxf database.tgz;
rm database.tgz;
#########################################
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/disulfidize_protein/run.sh;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/disulfidize_protein/flags;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/disulfidize_protein/disulfidize.xml;
# run disulfidize script
sh run.sh $disulfide.pdb &
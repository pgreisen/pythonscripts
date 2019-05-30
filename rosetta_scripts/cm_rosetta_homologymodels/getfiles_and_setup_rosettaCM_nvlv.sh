###############################################
hhsearchdir=$1;
aws s3 cp s3://tempfilescm/$hhsearchdir.zip . 
unzip $hhsearchdir.zip;
rm $hhsearchdir.zip;
###############################################
mkdir bin;
cd bin;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/cm_rosetta_homologymodels/setup_RosettaCM.py;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/cm_rosetta_homologymodels/get_top_n_sequence.py;
mkdir rosetta_bin;
cd rosetta_bin;
# test if files have already been copied
if [! -f "partial_thread.static.linuxgccrelease" ]; then
    aws s3 cp s3://enevolvcomputationalbiology/programs/partial_thread.static.linuxgccrelease.tgz .
fi
if [! -f "rosetta_scripts.static.linuxgccrelease" ]; then
    aws s3 cp s3://enevolvcomputationalbiology/programs/rosetta_scripts.static.linuxgccrelease.tgz .
fi
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
cd ../..;
if [ ! -d "database" ]; then
    aws s3 cp s3://enevolvcomputationalbiology/programs/database.tgz .;
    tar zxf database.tgz;
    rm database.tgz;
fi
#########################################
cd $hhsearchdir;
fasta=`find *fasta`;
hhsearch=`find *_pdb.hhr`;
python ~/bin/get_top_n_sequence.py $hhsearch 10;
python2.7 ~/bin/setup_RosettaCM.py --fasta $fasta --alignment reduced_hhsearch.hhr --alignment_format hhsearch --rosetta_bin ~/bin/rosetta_bin --build static;
############################################
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/cm_rosetta_homologymodels/qsub.sh;
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/cm_rosetta_homologymodels/run_setup.sh;
pth=$hhsearchdir/rosetta_cm
for i in $(seq 1 100);
do
    qsub -v path=$pth,prefix=$i qsub.sh;
done

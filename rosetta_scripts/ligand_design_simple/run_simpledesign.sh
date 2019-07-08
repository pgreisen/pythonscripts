#!/bin/bash -f 
# Get rosetta executables and database
pth_on_S3=s3://enevolvcomputationalbiology/programs;
pth_on_S3_pdbs=s3://temppdb/;
# update machine to contain new libraries
sudo apt update -y
DEBIAN_FRONTEND=noninteractive apt-get upgrade -yq
sudo apt install -y build-essential
sudo apt-get install -y libsqlite3-dev
sudo apt install -y libgl1-mesa-dev
sudo apt  install -y awscli 
sudo apt install -y emacs
sudo apt install -y unzip
sudo apt install -y zip
sudo apt-get install -y zlib1g-dev 
sudo apt-get install -y scons 
pth=`pwd`
if [ ! -d "database" ]; then
    aws s3 cp $pth_on_S3/database.tgz .;
    tar zxf database.tgz;
    database=$pth/database;
    rm database.tgz;
fi

if [ ! -f "rosetta_scripts.static.linuxgccrelease" ]; then
    ros=rosetta_scripts.static.linuxgccrelease;
    aws s3 cp $pth_on_S3/$ros.tgz . ;
    tar zxf $ros.tgz;
    rm $ros.tgz;
    exe=$pth/$ros;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_design_simple/flags;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ligand_design_simple/ligandesign_simple.xml;
fi

# copy files over
aws s3 cp $pth_on_S3_pdbs . --include="*.zip" --recursive;
# generalize it

# setup and run design
nthreads=5;
nstruct=20;
pth=`pwd`;
rosetta_db=$pth/database;
for k in *.fa.params;
do
    params=$k;
done

for i in *.pdb;
do
    for j in $(seq 1 $nthreads);
    do
	nice ~/rosetta_scripts.static.linuxgccrelease -database $rosetta_db -extra_res_fa $params -parser:protocol $pth/ligandesign_simple.xml -in:file:s $i -in:file:native $i -nstruct $nstruct -out:file:scorefile liganddesign.sc -packing:ex1 -packing:ex2 -out:prefix $j\_simple @$pth/flags > /dev/null & echo "done";
    done
done

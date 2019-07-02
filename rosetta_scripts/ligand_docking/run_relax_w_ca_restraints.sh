#!/bin/bash -f

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
    git clone https://github.com/pgreisen/pythonscripts/tree/master/rosetta_scripts/ligand_docking;
    mv ligand_docking/* .;
    rm -rf ligand_docking;
fi

# copy files over
aws s3 cp $pth_on_S3_pdbs . --include="*.zip" --recursive;
# generalize it
for tmp in *zip;
do
    unzip $tmp;
    rm $tmp;
done

# zip file must contain pdb files
for i in *pdb;
do
    dst=${i%.pdb};
    mkdir $dst;
    mv $i $dst/;
    cd $dst;
    run_relax & echo "Running job with pdb $i";
    cd ..;
done
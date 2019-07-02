#!/bin/bash -f

run_relax() {
    i=1;
    ncycles=10;
    ddgcut=1.0;
    start=3;
while [ $i -lt $ncycles ]
do
    newpdb=`ls -t *pdb | head -1`;
    echo $newpdb;
    $exe -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -s $newpdb @$pth/flags -database $database > /dev/null;
    initscore=$(awk 'NR == start' start=$start score.sc | awk '{print $2}');
    echo $initscore;
    if [ $i -gt 1 ]
    then
        num=$((start + i - 1));
        newscore=$(awk 'NR == num' num=$num score.sc | awk '{print $2}');
        if [ $i -gt 2 ]
        then
            num=$((start + i - 2));
            initscore=$(awk 'NR == num' num=$num score.sc | awk '{print $2}');
        fi
	echo "Newscore: $newscore";
        delta=$(echo "$initscore - $newscore" | bc);
        echo "Deltascore: $delta";
        evalend="$(echo "${delta} < ${ddgcut}" | bc)";
        if [ 1 -eq  "$evalend" ]
        then
            break;
        fi
    fi
    i=$(($i+1))
done
}

# Get rosetta executables and database
pth_on_S3=s3://enevolvcomputationalbiology/programs;
pth_on_S3_pdbs=s3://temppdb/;
pth=`pwd`
if [ ! -d "database" ]; then
    aws s3 cp $pth_on_S3/database.tgz .;
    tar zxf database.tgz;
    database=$pth/database;
fi

if [ ! -f "relax.static.linuxgccrelease" ]; then
    aws s3 cp $pth_on_S3/relax.static.linuxgccrelease.tgz . ;
    tar zxf relax.static.linuxgccrelease.tgz;
    exe=$pth/relax.static.linuxgccrelease;
    wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/relax_w_ca_restraints/flags;
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
    dst=${i%.pdb}
    mv $i $dst;
    cd $dst;
    run_relax;
    cd ..;
done
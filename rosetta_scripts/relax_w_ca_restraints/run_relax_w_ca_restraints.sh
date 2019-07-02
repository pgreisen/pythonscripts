#!/bin/bash -f
# Get rosetta executables and database
pth_on_S3=s3://enevolvcomputationalbiology/programs;
pth_on_S3_pdbs=s3://temppdb/;
if [ ! -d "database" ]; then
    aws s3 cp $pth_on_S3/database.tgz .;
    tar zxf database.tgz;
    database=/home/ubuntu/database;
fi

if [ ! -f "relax.static.linuxgccrelease" ]; then
    aws s3 cp $pth_on_S3/relax.static.linuxgccrelease.tgz . ;
    tar zxf relax.static.linuxgccrelease.tgz;
    exe=/home/ubuntu/relax.static.linuxgccrelease;
fi

# copy files over
aws s3 cp $pth_on_S3_pdbs . --include="*.zip" --recursive;

for tmp in *zip;
do
    unzip $tmp;
    rm $tmp;
done


i=1;
ncycles=10;
ddgcut=1.0;

start=3;
while [ $i -lt $ncycles ]
do
    newpdb=`ls -t *pdb | head -1`;
    echo $newpdb;
    $exe -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -s $newpdb -database $database > /dev/null;
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

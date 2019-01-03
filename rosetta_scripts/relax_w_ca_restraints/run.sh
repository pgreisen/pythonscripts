#!/bin/bash -f
i=1;
ncycles=8;
ddgcut=1.0;
exe=/home/ubuntu/relax.static.linuxgccrelease;
database=/home/ubuntu/database;
# cd $PBS_O_WORKDIR;
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
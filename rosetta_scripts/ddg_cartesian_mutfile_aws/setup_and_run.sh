# First minimize the initial pdb
# it is required for the scripts to work
# that rosetta database and executables
# are linked or exists from home
#
######################
##     Enevolv
######################
aws s3 cp s3://enevolvcomputationalbiology/programs/cartesian_ddg.static.linuxgccrelease.tgz .
aws s3 cp s3://enevolvcomputationalbiology/programs/relax.static.linuxgccrelease.tgz .
aws s3 cp s3://enevolvcomputationalbiology/programs/database . --recursive --exclude "*" --include "*"
####################
initpdb=$1;
mindir=min_cart;
mkdir $mindir;
cd $mindir;
echo "switch:cartesian
repeat 2
ramp_repack_min 0.02  0.01     1.0  50
ramp_repack_min 0.250 0.01     0.5  50
ramp_repack_min 0.550 0.01     0.0 100
ramp_repack_min 1     0.00001  0.0 200
accept_to_best
endrepeat" > cart2.script;
ln -s ../$initpdb 1.pdb;
sh ../min_cart.sh;
pdbfile=`sort -n -k2 score.sc | awk '{print $NF}' | head -n 1`;
cp $pdbfile.pdb ../lowest.pdb;
cd ..;
# We loop over all resfile and setup directories for 
# their calculations
for i in *.mutfile;
do
    dir=${i%*.mutfile};
    mkdir $dir;
    mv $i $dir/;
    cd $dir;
    for j in mutfile_*;
    do
        resf=$j;
    done

    pth=`pwd`;
    echo $pth;
    nohup qsub -v pth=$pth,mutfile=$resf ~/qsub_run_ddg.sh & echo "done";
    cd ..;
done

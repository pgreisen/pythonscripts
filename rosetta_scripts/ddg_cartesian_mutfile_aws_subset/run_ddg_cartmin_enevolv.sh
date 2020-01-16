##############################################




######################################
# copy script to update ec2 instance - assume it is ubuntu for now
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws_subset/update_ec2_ubuntu_aws_no_python.sh
sh update_ec2_ubuntu_aws_no_python.sh;
####################################### First minimize the initial pdb
# it is required for the scripts to work
# that rosetta database and executables
# are linked or exists from home
#
######################
##     set path to S3 bucket
######################
pth_on_S3=s3://enevolvcomputationalbiology/programs;
pth_on_S3_pdbs=s3://tmppdb/;


aws s3 cp $pth_on_S3/cartesian_ddg.static.linuxgccrelease.tgz .
aws s3 cp $pth_on_S3/relax.static.linuxgccrelease.tgz .
aws s3 cp $pth_on_S3/database.tgz .;

# unpack files
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
# Next, get pdb file on S3
aws s3 cp $pth_on_S3_pdbs . --include="*.zip" --recursive;

for tmpzip in *zip;
do 
    unzip $tmpzip;
done

####################
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ddg_cartesian_mutfile_aws_subset/qsub_run_ddg.sh

for pdb in *.pdb;
do
    tmppdb=$pdb;
done

cp $tmppdb 1.pdb;
echo $tmppdb;
initpdb=1.pdb;
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
mv ../1.pdb .;
# start the cartesian minimization before substitution
nohup ~/relax.static.linuxgccrelease -s 1.pdb -use_input_sc -constrain_relax_to_start_coords -ignore_unrecognized_res -nstruct 20 -relax:coord_constrain_sidechains -relax:cartesian -score:weights ref2015_cart -relax:min_type lbfgs_armijo_nonmonotone -relax:script cart2.script -database ~/database;
#####
pdbfile=`sort -n -k2 score.sc | awk '{print $NF}' | head -n 1`;
cp $pdbfile.pdb ~/lowest.pdb;
cd ..;
# We loop over all resfile and setup directories for 
# their calculations
for i in *.mutfile;
do
    dir=${i%*.mutfile};
    mkdir $dir;
    mv $i $dir/;
    cd $dir;
    cp ../lowest.pdb .;
    for j in *.mutfile;
    do
        resf=$j;
    done
    pth=`pwd`;
    nohup qsub -v pth=$pth,mutfile=$resf ~/qsub_run_ddg.sh & echo "done";
    cd ..;
done

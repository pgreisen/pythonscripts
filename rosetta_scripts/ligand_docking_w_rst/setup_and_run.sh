# First minimize the initial pdb
# it is required for the scripts to work
# that rosetta database and executables
# are linked or exists from home
#
######################
##     Shenzhen
######################
aws s3 cp s3://programexecutables/rosetta_scripts.static.linuxgccrelease.tgz .
aws s3 cp s3://programexecutables/database.tgz .;
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
####################
initpdb=$1;
ln -s $initpdb 1.pdb;
params=`find *params`;
nprc=`nproc`;

for i in $(seq 1 $nprc);
do
    sh run_ligand_docking.sh $params $initpdb $i > /dev/null & echo "Done";
done

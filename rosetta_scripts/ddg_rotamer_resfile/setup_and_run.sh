# First minimize the initial pdb
# it is required for the scripts to work
# that rosetta database and executables
# are linked or exists from home
#
initpdb=$1;
mkdir min_w_cst;
cd min_w_cst;
echo 1.pdb > lst;
ln -s ../$initpdb 1.pdb;
sh ../run_min_cst.sh;
cd ..;
# We loop over all resfile and setup directories for 
# their calculations
for i in resfile*;
do
    dir=${i#*_};
    mkdir $dir;
    mv $i $dir/;
    cd $dir;
    for j in resfile*;
    do
        resf=$j;
    done
    ln -s $resf resfile;
    nohup sh ../run_ddg.sh & echo "done";
    cd ..;
done
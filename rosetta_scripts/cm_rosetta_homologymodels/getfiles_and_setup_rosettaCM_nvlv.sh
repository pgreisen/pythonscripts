mkdir bin;
mv setup_RosettaCM.py bin/;
mv get_top_n_sequence.py bin/;
cd bin;
mkdir rosetta_bin;
cd rosetta_bin;
aws s3 cp s3://programexecutables/partial_thread.static.linuxgccrelease.tgz .
aws s3 cp s3://programexecutables/rosetta_scripts.static.linuxgccrelease.tgz .
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done
cd ../..;
aws s3 cp s3://programexecutables/database.tgz .;
tar zxf database.tgz;
rm database.tgz;
cd ;
cp qsub.sh *rosetta_cm/;
sh run_setup.sh & echo "Running";
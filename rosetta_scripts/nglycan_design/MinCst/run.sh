# run program
exe=~/rosetta_bin/rosetta_scripts.static.linuxgccrelease;
database=~/database;
xml=./min_cst.xml;
flags=./flags;
# run program
for i in *pdb;
do
    $exe -database $database -parser:protocol $xml -in:file:s $i @$flags
done
exit
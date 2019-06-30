# run program
exe=~/rosetta_scripts.static.linuxgccrelease;
database=~/database;
xml=./min_cst.xml;
flags=./flags;
# run program
for i in *;
do
    $exe -database $database -parser:protocol $xml -in:file:s $i @$flags > /dev/null;
done
exit

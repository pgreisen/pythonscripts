# run program
exe=~/rosetta_bin/rosetta_scripts.static.linuxgccrelease;
database=./database;
xml=./docking_minimize.xml;
flags=./flags;
# run program
$exe -database $database -parser:protocol $xml -in:file:s $1 @$flags

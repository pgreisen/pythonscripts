# run program
exe=../ScanDisulfides/rosetta_scripts.static.linuxgccrelease;
database=~/database;
xml=./hbnet.xml;
flags=./flags;
# run program
$exe -database $database -parser:protocol $xml -in:file:s $1 @$flags
exit

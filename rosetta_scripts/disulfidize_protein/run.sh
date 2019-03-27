# run program
exe=~/rosetta_bin/rosetta_scripts.static.linuxgccrelease;
database=~/database;
xml=~/disulfidize.xml;
flags=~/flags;
# run program
$exe -database $database -parser:protocol $xml -in:file:s $1 @$flags
exit

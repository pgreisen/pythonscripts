# run program
exe=./rosetta_scripts.static.linuxgccrelease;
database=~/database;
xml=disulfidize_across_interface.xml
flags=./flags;
# run program
$exe -database $database -parser:protocol $xml -in:file:s $1 @$flags
exit

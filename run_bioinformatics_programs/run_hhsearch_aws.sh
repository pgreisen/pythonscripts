export HHLIB=/home/ubuntu/hh-suite;
exe_hhsearch=/home/ubuntu/hh-suite/build/bin/hhsearch;
exe_hhblist=/home/ubuntu/hh-suite/build/bin/hhblits;
database_seqs=/home/ubuntu/hhsearch_databases/pfam/pfam;
database_pdb=/home/ubuntu/hhsearch_databases/pdb70;
####################
# gather seequences
$exe_hhblist -d $database_seqs -i $1 -v 1 -oa3m ${1%.fasta}.a3m -n 3
# gather homologous structures
$exe_hhsearch -i ${1%.fasta}.a3m -d $database_pdb -o ${1%.fasta}\_pdb.hhr

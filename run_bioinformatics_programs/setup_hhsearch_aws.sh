export HHLIB=/home/ubuntu/hh-suite;
exe_hhsearch=/home/ubuntu/build/bin/hhsearch;
exe_hhblist=/home/ubuntu/build/bin/hhblits;
database_seqs=/home/ubuntu/hhsearch_databases/
# setup database
cd $database_seqs;
mkdir pfam;
mv pfamA_32.0.tar.gz pfam;
cd pfam;
tar zxf pfamA_32.0.tar.gz & cd ..;
mkdir pdb70;
mv pdb70_from_mmcif_latest.tar.gz pdb70/;
cd pdb70;
tar zxf pdb70_from_mmcif_latest.tar.gz & cd ..;



$exe_hhblist -d /databases/pfam/pfam/pfam -i $1 -v 1 -oa3m ${1%.fasta}.a3m -n 3
#
/tools/hh-suite/build/bin/hhsearch -i ${1%.fasta}.a3m -d /databases/pdb70/pdb70 -o ${1%.fasta}\_pdb.hhr

sudo docker cp 9659b2aa6ee3:/root/test.fasta foo.txt
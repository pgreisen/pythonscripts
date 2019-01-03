######################################
# copy script to update ec2 instance - assume it is ubuntu for now
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws/update_ec2_ubuntu_aws_no_python.sh
sh update_ec2_ubuntu_aws_no_python.sh;
######################################
# next copy executables from S3 bucket
aws s3 cp s3://enevolvcomputationalbiology/databases/hh-suite.tgz .
tar zxf hh-suite.tgz;
rm hh-suite.tgz;
aws s3 cp s3://enevolvcomputationalbiology/databases/hhsearch_databases.tgz .
tar zxf hhsearch_databases.tgz 
rm hhsearch_databases.tgz;
######################################
# get fasta files to work on
aws s3 cp s3://tempfilespssm . --include="*.zip" --recursive
######################################
export HHLIB=/home/ubuntu/hh-suite;
exe_hhsearch=/home/ubuntu/hh-suite/build/bin/hhsearch;
exe_hhblist=/home/ubuntu/hh-suite/build/bin/hhblits;
database_seqs=/home/ubuntu/hhsearch_databases/pfam/pfam;
database_pdb=/home/ubuntu/hhsearch_databases/pdb70;
####################
# gather seequences
unzip *zip;
rm *zip;
for i in *.fasta;
do
    dst=${i%.fasta}\_hhsearch;
    mkdir $dst;
    mv $i $dst;
    cd $dst;
    $exe_hhblist -d $database_seqs -i $1 -v 1 -oa3m ${1%.fasta}.a3m -n 3
    # gather homologous structures
    $exe_hhsearch -i ${1%.fasta}.a3m -d $database_pdb -o ${1%.fasta}\_pdb.hhr
    cd ..;
done




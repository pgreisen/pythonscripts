######################################
# copy script to update ec2 instance - assume it is ubuntu for now
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws/update_ec2_ubuntu_aws_no_python.sh
sh update_ec2_ubuntu_aws_no_python.sh;
######################################
# next copy executables from S3 bucket
# test if files have already been copied
if [ ! -d "hh-suite" ]; then
    aws s3 cp s3://proteindatabases//hh-suite.tgz .;
    tar zxf hh-suite.tgz;
    rm hh-suite.tgz;
fi
if [ ! -d "hhsearch_databases" ]; then
    aws s3 cp s3://proteindatabases/hhsearch_databases_2019.tgz .
    tar zxf hhsearch_databases_2019.tgz 
    rm hhsearch_databases_2019.tgz;
    # unpack files
    cd hhsearch_databases;
    mkdir pfam;
    mv pfamA_32.0.tar.gz pfam;
    cd pfam;
    tar zxf pfamA_32.0.tar.gz;
    rm pfamA_32.0.tar.gz;
    cd ..;
    mkdir pdb70;
    mv pdb70_from_mmcif_latest_2019.tar.gz pdb70/;
    cd pdb70;
    tar zxf pdb70_from_mmcif_latest_2019.tar.gz; 
    rm pdb70_from_mmcif_latest_2019.tar.gz; 
    cd ;
fi
######################################
######################################
# get fasta files to work on
aws s3 cp s3://tmpfilefasta . --include="*.zip" --recursive;
######################################
pth=`pwd`

export HHLIB=$pth/hh-suite;
exe_hhsearch=$pth/hh-suite/build/bin/hhsearch;
exe_hhblist=$pth/hh-suite/build/bin/hhblits;
database_seqs=$pth/hhsearch_databases/pfam/pfam;
database_pdb=$pth/hhsearch_databases/pdb70/pdb70;
####################
# gather seequences from fasta file
unzip *zip;
rm *zip;
for i in *.fasta;
do
    dst=${i%.fasta}\_hhsearch;
    mkdir $dst;
    mv $i $dst;
    cd $dst;
    $exe_hhblist -d $database_seqs -i $i -v 1 -oa3m ${i%.fasta}.a3m -n 3
    # gather homologous structures
    $exe_hhsearch -i ${i%.fasta}.a3m -d $database_pdb -o ${i%.fasta}\_pdb.hhr
    cd ..;
done
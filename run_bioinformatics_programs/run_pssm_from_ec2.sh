######################################
# copy script to update ec2 instance - assume it is ubuntu for now
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws/update_ec2_ubuntu_aws_no_python.sh
sh update_ec2_ubuntu_aws_no_python.sh;
######################################
# next copy executables from S3 bucket
aws s3 cp s3://tmpfilefasta . --include="*.zip" --recursive;
unzip *zip;
rm *zip;
aws s3 cp s3://proteindatabases/uniref90files.zip .;
######################################
# setup for pssm run
path=~/ncbi-blast-2.7.1+/bin
database=~/uniref90
######################################
# unzip file for pssm
unzip uniref90files.zip;
rm uniref90files.zip;
for i in *.fasta;
do
    dst=${i%.fasta}\_pssm;
    mkdir $dst;
    mv $i $dst;
    cd $dst;
    $path/psiblast -query $i -db $database/uniref90_w_index.fasta -out_pssm my_protein.ckp -evalue 0.01 -out_ascii_pssm ascii_mtx_file -out output_file -num_iterations 3 > pssm_logfile & echo "done";
    cd ..;
done
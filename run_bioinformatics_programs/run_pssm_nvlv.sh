######################################
# copy script to update ec2 instance - assume it is ubuntu for now
wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws/update_ec2.sh && bash update_ec2.sh;
######################################
# next copy executables from S3 bucket
aws s3 cp s3://tempfilespssm . --include="*.zip" --recursive;
for i in *zip;
do
    unzip $i;
    rm $i;
done
# setup the database necessary for the run
pth=`pwd`;
path=$pth/ncbi-blast-2.7.1+/bin
database=$pth/uniref90

if [ ! -d "uniref90" ]; then
    aws s3 cp s3://enevolvcomputationalbiology/databases/uniref90files.zip .;
    unzip uniref90files.zip;
    rm uniref90files.zip;
fi

for i in *.fasta;
do
    dst=${i%.fasta}\_pssm;
    mkdir $dst;
    mv $i $dst;
    cd $dst;
    $path/psiblast -query $i -db $database/uniref90_w_index.fasta -out_pssm my_protein.ckp -evalue 0.01 -out_ascii_pssm ascii_mtx_file -out output_file -num_iterations 3 > pssm_logfile & echo "done";
    cd ..;
done

######################
##     set path to S3 bucket
######################
pth_on_S3=s3://programexecutables;
pth_on_S3_pdbs=s3://tmppdb/;

wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/bashscripts/update_ec2_aws/update_ec2_ubuntu_aws_no_python.sh;



aws s3 cp $pth_on_S3/loopmodel.static.linuxgccrelease.tgz .
aws s3 cp $pth_on_S3/database.tgz .;
# unpack files
for i in *tgz;
do
    tar zxf $i;
    rm $i;
done


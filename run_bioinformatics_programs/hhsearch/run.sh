# input is fasta file and this is setup for my account as an initial test
echo "Initial input is the fastafile - make sure the path to the credentials has been set correctly";
runfile=/Users/pgreisen/Projects/Enevolv/code/enevolvcode/BioinformaticAnalysis/HHsearch/run_hhsearch_aws_nvlv.sh;
keyname='/Users/pgreisen/Downloads/compilerosetta.pem';
fastafile=test.fasta;
aws_credentials=/Users/pgreisen/Projects/Enevolv/code/enevolvcode/AWS/credentials/credentials.txt
zip ${fastafile%.fasta}.zip $fastafile;
fastazipped=${fastafile%.fasta}.zip;
bucket=tempfilespssm;
pth=/Users/pgreisen/Projects/Enevolv/code/enevolvcode/BioinformaticAnalysis/HHsearch;
awsupdates=/Users/pgreisen/Projects/Enevolv/code/enevolvcode/AWS/update_ec2_ubuntu/update_ec2.sh;

python $pth/aws_spinup_and_execute_hhsearch.py -z $fastazipped -f $aws_credentials -b $bucket -r $runfile -k $keyname --awsupdate $awsupdates

# input is fasta file and this is setup for my account as an initial test
echo "Initial input is the fastafile - make sure the path to the credentials has been set correctly";
runfile='https://raw.githubusercontent.com/pgreisen/pythonscripts/master/run_bioinformatics_programs/run_pssm_from_ec2.sh';
keyname='/Users/pgreisen/AWS/tmp_keys/4TestKey.pem';
fastafile=SLC16A11.fasta;
zip ${fastafile%.fasta}.zip $fastafile;
fastazipped=${fastafile%.fasta}.zip;

python ../aws_spinup_and_execute_pssm.py -z $fastazipped -f ~/.aws/credentials -b tmpfilefasta -r $runfile -k $keyname

# input is fasta file and this is setup for my account as an initial test
echo "Initial input name of cluster and configuration file for setting up cluster";
configfile=/Users/pgreisen/.cfncluster/config;
keyname='/Users/pgreisen/AWS/tmp_keys/4TestKey.pem'
aws_credentials=/Users/pgreisen/AWS/credentials/credentials;
ec2update=/Users/pgreisen/pythonscripts/bashscripts/update_ec2_aws/update_ec2.sh;
python aws_spinup_and_create_cluster.py -f $aws_credentials -k $keyname --cluster_config_file $configfile --ec2_update $ec2update

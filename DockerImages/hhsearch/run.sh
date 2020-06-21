runfile=/home/run_hhsearch_aws_nvlv.sh
keyname=`find /home/*.pem`
fastazipped=/home/fastafile.zip
bucket=tempfilespssm
python /home/aws_spinup_and_execute_hhsearch.py -z $fastazipped -b $bucket -r $runfile -k $keyname --awsupdate $awsupdates /home/update_ec2.sh --aws_secret_access_key $AWS_SECRET_ACCESS_KEY --aws_access_key_id $AWS_ACCESS_KEY_ID

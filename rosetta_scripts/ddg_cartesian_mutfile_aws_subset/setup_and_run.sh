###########################################################################
#
#
# get the files from S3 and setup the runs
# Make sure that the credentials have been set correctly
#
#
###########################################################################
nohup wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/rosetta_scripts/ddg_cartesian_mutfile_aws_subset/run_ddg_cartmin.sh && sh run_ddg_cartmin.sh;

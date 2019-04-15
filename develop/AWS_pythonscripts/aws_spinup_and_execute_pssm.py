import boto3
import paramiko,time
import argparse
from boto3.s3.transfer import S3Transfer

class RunPSSM:

    def __init__(self):
        self.instance_type = 't2.xlarge'
        self.credentials = ""
        self.file_w_credentials = ""
        self.zipped_fasta_file = ""
        self.imageid = 'ami-0ac019f4fcb7cb7e6'
        # Gb on running instance
        self.volumensize = 150
        self.keyname = ""
        self.file2copy = ""
        self.bucket = ""
        self.aws_access_key_id = ""
        self.aws_secret_access_key = ""


    def ssh_and_execute_script(self,host,username,command):
        # time for aws to spin up the instance
        time.sleep(30)
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(hostname=host,username=username,timeout=60)
        print("Connected to the server",host)
        print("Executing command --> {}".format(command))
        stdin, stdout, stderr = client.exec_command(command,timeout=10)
        ssh_output = stdout.read()
        ssh_error = stderr.read()
        if ssh_error:
            print("Problem occurred while running command:"+ command + " The error is " + ssh_error)
        else:
            print("Job running")

    def create_ec2_instance(self,):
        ec2 = boto3.resource('ec2')
        # create a new EC2 instance
        instances = ec2.create_instances(
            ImageId=self.imageid,
            MinCount=1,
            MaxCount=1,
            InstanceType=self.instance_type,
            KeyName=self.key,
            BlockDeviceMappings=[{"DeviceName": "/dev/sda1","Ebs" : { "VolumeSize" : self.volumensize }}])

    def copy_file_to_S3(self):

        session = boto3.Session(aws_access_key_id=self.aws_access_key_id, aws_secret_access_key=self.aws_secret_access_key)
        s3 = session.resource('s3')
        client = s3.meta.client
        transfer = S3Transfer(client)
        transfer.upload_file(self.zipped_fasta_file, self.bucket, self.zipped_fasta_file)


    def set_credentials(self):
        with open(self.credentials,'r') as f:
            for line in f:
                if("aws_access_key_id" in line):
                    self.aws_access_key_id = self.aws_access_key_id+line.strip().split("=")[1].strip()
                elif("aws_secret_access_key" in line):
                    self.aws_secret_access_key += line.strip().split("=")[1].strip()
                else:
                    continue


    def main(self):

        parser = argparse.ArgumentParser(description="Setup PSSM calculation using AWS")
        parser.add_argument('-i', dest='instancetype', help='Compute instance to use for calculation')
        # parser.add_argument('-f', dest='file_w_credentials', help='File containing the credentials used to transfer files between S3 and EC2')
        parser.add_argument('-z', dest='zipped_fasta_file', help='Zip-file to transfer')
        parser.add_argument('-m', dest='imageid', help='ImageId', default="ami-0ac019f4fcb7cb7e6")
        parser.add_argument('-c', dest='credentials', help='')
        parser.add_argument('-b', dest='bucket', help='')
        parser.add_argument('-k', dest='keyname', help='')




        args_dict = vars(parser.parse_args())

        for item in args_dict:
            setattr(self, item, args_dict[item])
        # getting credential
        self.set_credentials()
        #
        self.copy_file_to_S3()




if __name__ == "__main__":
    run = RunPSSM()
    run.main()







'''
# Boto 3
# Use the filter() method of the instances collection to retrieve
# all running EC2 instances.
running_instances = ec2.instances.filter(
    Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
for instance in running_instances:
    print(instance.id, instance.instance_type)

url_login=instance.public_dns_name
ssh_and_execute_script(url_login,"ubuntu","wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/run_bioinformatics_programs/run_pssm_nvlv.sh && sh run_pssm_nvlv.sh")
'''
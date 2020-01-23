import boto3
import paramiko,time
import argparse
from boto3.s3.transfer import S3Transfer
import shutil, os

class RunPSSM:

    def __init__(self):
        self.instance_type = 't2.xlarge'
        self.username = "ubuntu"
        self.credentials = ""
        self.file_w_credentials = ""
        self.zipped_fasta_file = ""
        # a public image id for ubuntu
        self.imageid = 'ami-0ac019f4fcb7cb7e6'
        # Gb on running instance
        self.volumensize = 150
        self.key = ""
        self.file2copy = ""
        self.bucket = ""
        self.aws_access_key_id = ""
        self.aws_secret_access_key = ""
        self.session = ""
        self.instance_id = ""
        self.url_instance = ""
        self.keyname = ""
        self.run_file = ""
        self.exe_pssm = ""
        # time to wait for the ec2 instance to spin up
        self.sleep = 60
        self.vpc = "vpc-0a3c8ad289df9821e"
        self.subnet = "subnet-03c72665ab3e78032"
        self.sg = "sg-0a88bfaa0bc258091"

    def set_run_file(self):
        file = self.run_file.split('/')[-1]
        #self.exe_pssm = "nohup bash "+file+" /dev/null 2>&1 &"
        self.exe_pssm = str("nohup bash " + file + " & echo ")


    def get_setup_and_run(self):
        getfile = "wget "+self.run_file
        self.set_run_file()
        exe='''
        mkdir .aws;
        echo [default] >>.aws/credentials \n echo aws_access_key_id = '''+self.aws_access_key_id+''' >>.aws/credentials \n echo aws_secret_access_key = '''+self.aws_secret_access_key+\
            ''' >> .aws/credentials \n'''+getfile+''' && '''+self.exe_pssm+''';'''
        print("Type of return from get_setup_and_run: ", type(exe), exe)
        return exe

    def ssh_and_execute_script(self, host, username, command):
        print(self.keyname)
        cert = paramiko.RSAKey.from_private_key_file(self.keyname)

        # time for aws to spin up the instance
        time.sleep(self.sleep)
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(hostname=host, username=username, timeout=self.sleep, pkey = cert)
        print("Connected to the server",host)
        print("Executing command --> {}".format(command))
        stdin, stdout, stderr = client.exec_command(command, timeout=self.sleep)
        ssh_output = stdout.read()
        print(ssh_output)
        ssh_error = stderr.read()
        if ssh_error:
            print("#####1: ",command)
            print("#########2: ", ssh_error)
            print("Problem occurred while running command:"+ command + " The error is " + ssh_error)
        else:
            print("Job running")

    def create_ec2_instance(self):
        print(self.vpc)
        self.set_credentials()
        ec2 = self.session.resource('ec2')
        vpc = ec2.Vpc(self.vpc)
        # Attributes
        print(vpc.cidr_block)
        print(vpc.state)

        # associate the route table with the subnet
        self.set_keyname()
        print(self.key,ec2,type(self.key))
        # create a new EC2 instance

        instances = ec2.create_instances(
            ImageId=self.imageid,
            MinCount=1,
            MaxCount=1,
            InstanceType=self.instance_type,
            KeyName=str(self.key),
            BlockDeviceMappings=[{"DeviceName": "/dev/sda1","Ebs" : { "VolumeSize" : self.volumensize }}],
            NetworkInterfaces=[{'SubnetId': self.subnet, 'DeviceIndex': 0, 'AssociatePublicIpAddress': True,'Groups': [self.sg]}])

        for i in instances:
            self.instance_id = i.id
        #
        instances[0].wait_until_running()
        print(instances[0].id)

        from time import sleep
        # time to setup ec2 instance
        sleep(self.sleep)
        running_instances = ec2.instances.filter(Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
        for instance in running_instances:
            print(self.instance_id,instance.id)
            if(self.instance_id == instance.id ):
                self.url_instance = instance.public_dns_name
        print(self.url_instance)


    def copy_file_to_S3(self) -> object:
        session = boto3.Session(aws_access_key_id=self.aws_access_key_id, aws_secret_access_key=self.aws_secret_access_key)
        s3 = session.resource('s3')
        client = s3.meta.client
        transfer = S3Transfer(client)
        transfer.upload_file(self.zipped_fasta_file, self.bucket, self.zipped_fasta_file)

    def set_keyname(self):
        self.key = self.keyname.split("/")[-1].split(".")[0]

    def set_credentials_constructor(self):
        with open(self.file_w_credentials ,'r') as f:
            for line in f:
                if("aws_access_key_id" in line):
                    self.aws_access_key_id = self.aws_access_key_id+line.strip().split("=")[1].strip()
                elif("aws_secret_access_key" in line):
                    self.aws_secret_access_key += line.strip().split("=")[1].strip()
                else:
                    continue

    def set_credentials(self):
        self.session = boto3.Session(aws_access_key_id=self.aws_access_key_id,aws_secret_access_key=self.aws_secret_access_key)

    def main(self):

        parser = argparse.ArgumentParser(description="Setup PSSM calculation using AWS")
        parser.add_argument('-i', dest='instancetype', help='Compute instance to use for calculation')
        parser.add_argument('-f', dest='file_w_credentials', help='File containing the credentials used to transfer files between S3 and EC2')
        parser.add_argument('-z', dest='zipped_fasta_file', help='Zip-file to transfer')
        parser.add_argument('-m', dest='imageid', help='ImageId', default="ami-0ac019f4fcb7cb7e6")
        parser.add_argument('-c', dest='credentials', help='')
        parser.add_argument('-b', dest='bucket', help='Destination of the zip file used to run the desired application')
        parser.add_argument('-k', dest='keyname', help='Key used to access the instance')
        parser.add_argument('-r', dest='run_file', help='Bash script containing the run specifications')
        parser.add_argument('--vpc', dest='vpc', help='VPC - virtual private cloud to use for spinning up ec2 instances',default='vpc-0a3c8ad289df9821e')

        args_dict = vars(parser.parse_args())

        for item in args_dict:
            setattr(self, item, args_dict[item])

         # getting credential
        self.set_credentials_constructor()
        #
        self.copy_file_to_S3()
        # create ec2 instance and get url
        self.create_ec2_instance()

        command = self.get_setup_and_run()
        print(command, self.url_instance, self.username)
        self.ssh_and_execute_script(self.url_instance,self.username,command)


if __name__ == "__main__":
    run = RunPSSM()
    run.main()
import boto3
import paramiko,time
import argparse
from boto3.s3.transfer import S3Transfer
import shutil, os

class RunPSSM:

    def __init__(self):
        self.instance_type = 't2.small'
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
        # time to wait for the ec2 instance to spin up
        self.sleep = 60
        self.vpc = "vpc-0a3c8ad289df9821e"
        self.subnet = "subnet-03c72665ab3e78032"
        self.sg = "sg-0a88bfaa0bc258091"
        self.cfnclusterfile = ""
        self.clustername = ""
        self.cluster_settings = ""
        self.update_ec2_cmd = ""


    def get_update_ec2(self):
        updateec2 = ""
        with open(self.update_ec2_cmd, 'r') as f:
            for line in f:
                updateec2 += line
        return updateec2


    def set_clusterfile(self):
        with open(self.cfnclusterfile) as f:
            for line in f:
                self.cluster_settings += line


    def get_setup_and_run(self):
        self.set_clusterfile()
        exe='''
        bash update.sh;\nwait;\nmkdir .aws;
        echo [default] >>.aws/credentials \n echo aws_access_key_id = '''+self.aws_access_key_id+''' >>.aws/credentials \n echo aws_secret_access_key = '''+self.aws_secret_access_key+\
            ''' >> .aws/credentials \n mkdir .cfncluster; aws s3 cp s3://awsdependentfiles/cfncluster/config ~/.cfncluster/;sudo pip3 install cfncluster;cfncluster create '''+self.clustername
        return exe

    def ftp_files(self, client):
        updateec2 = self.get_update_ec2()
        ftp = client.open_sftp()
        file = ftp.file('update.sh', "a", -1)
        file.write(updateec2)
        file.flush()
        ftp.close()

    def ssh_and_execute_script(self, host, username, command):
        cert = paramiko.RSAKey.from_private_key_file(self.keyname)
        # time for aws to spin up the instance
        time.sleep(self.sleep)
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(hostname=host, username=username, timeout=self.sleep, pkey = cert)
        # Do ftp here?
        self.ftp_files(client)
        print("Connected to the server",host)
        print("Executing command --> {}".format(command))
        stdin, stdout, stderr = client.exec_command(command, timeout=self.sleep)
        ssh_output = stdout.read()
        ssh_error = stderr.read()
        if ssh_error:
            print("Problem occurred while running command:"+ command + " The error is " + ssh_error)
        else:
            print("Job running")
        # does this work?
        client.close()

    def create_ec2_instance(self):
        self.set_credentials()
        ec2 = self.session.resource('ec2')
        vpc = ec2.Vpc(self.vpc)
        # Attributes
        # associate the route table with the subnet
        self.set_keyname()
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

        from time import sleep
        # time to setup ec2 instance
        sleep(self.sleep)
        running_instances = ec2.instances.filter(Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
        for instance in running_instances:
            if(self.instance_id == instance.id ):
                self.url_instance = instance.public_dns_name

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

        parser = argparse.ArgumentParser(description="Create cfncluster on AWS")
        parser.add_argument('-i', dest='instancetype', help='Compute instance to use for calculation')
        parser.add_argument('-f', dest='file_w_credentials', help='File containing the credentials used to transfer files between S3 and EC2')
        parser.add_argument('-z', dest='zipped_fasta_file', help='Zip-file to transfer')
        parser.add_argument('-m', dest='imageid', help='ImageId', default="ami-0ac019f4fcb7cb7e6")
        parser.add_argument('-c', dest='credentials', help='')
        parser.add_argument('-b', dest='bucket', help='Destination of the zip file used to run the desired application')
        parser.add_argument('-k', dest='keyname', help='Key used to access the instance')
        parser.add_argument('-r', dest='run_file', help='Bash script containing the run specifications')
        parser.add_argument('-n','--name', dest='clustername', help='Name of cluster', default='cluster1')
        parser.add_argument('--cluster_config_file', dest='cfnclusterfile', help='file containing cluster settings')
        parser.add_argument('--ec2_update', dest='update_ec2_cmd', help='file containing how to update the ec2 instance')
        args_dict = vars(parser.parse_args())

        for item in args_dict:
            setattr(self, item, args_dict[item])

         # getting credential
        self.set_credentials_constructor()
        # create ec2 instance and get url
        self.create_ec2_instance()

        command = self.get_setup_and_run()
        self.ssh_and_execute_script(self.url_instance,self.username,command)


if __name__ == "__main__":
    run = RunPSSM()
    run.main()

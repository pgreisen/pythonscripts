import boto3
import paramiko
class spin_and_execute:
    
    def __init__(self):
        self.ssh_output = None
        self.ssh_error = None
        self.client = None
        #self.host= conf_file.HOST
        #self.username = conf_file.USERNAME
        #self.password = conf_file.PASSWORD
        #self.timeout = float(conf_file.TIMEOUT)
        #self.commands = conf_file.COMMANDS
        #self.pkey = conf_file.PKEY
        #self.port = conf_file.PORT
        #self.uploadremotefilepath = conf_file.UPLOADREMOTEFILEPATH
        #self.uploadlocalfilepath = conf_file.UPLOADLOCALFILEPATH
        #self.downloadremotefilepath = conf_file.DOWNLOADREMOTEFILEPATH
        #self.downloadlocalfilepath = conf_file.DOWNLOADLOCALFILEPATH


ec2 = boto3.resource('ec2')
# image ubuntu
## ami-0ac019f4fcb7cb7e6
# create a new EC2 instance
instances = ec2.create_instances(
     ImageId='ami-0ac019f4fcb7cb7e6',
     MinCount=1,
     MaxCount=1,
     InstanceType='t2.xlarge',
     KeyName='ec2-keypair_test2_boto3',
     BlockDeviceMappings=[{"DeviceName": "/dev/sda1","Ebs" : { "VolumeSize" : 50 }}]
)
# Boto 3
# Use the filter() method of the instances collection to retrieve
# all running EC2 instances.
running_instances = ec2.instances.filter(
    Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
for instance in running_instances:
    print(instance.id, instance.instance_type)

url_login=instance.public_dns_name

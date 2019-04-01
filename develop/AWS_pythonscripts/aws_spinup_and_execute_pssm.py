import boto3
import paramiko,time

def ssh_and_execute_script(host,username,command):
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
     BlockDeviceMappings=[{"DeviceName": "/dev/sda1","Ebs" : { "VolumeSize" : 180 }}]
)
# Boto 3
# Use the filter() method of the instances collection to retrieve
# all running EC2 instances.
running_instances = ec2.instances.filter(
    Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
for instance in running_instances:
    print(instance.id, instance.instance_type)

url_login=instance.public_dns_name
ssh_and_execute_script(url_login,"ubuntu","wget https://raw.githubusercontent.com/pgreisen/pythonscripts/master/run_bioinformatics_programs/run_pssm_nvlv.sh && sh run_pssm_nvlv.sh")

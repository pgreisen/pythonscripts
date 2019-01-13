import boto3

#chmod 400 ec2-keypair.pem

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
                                

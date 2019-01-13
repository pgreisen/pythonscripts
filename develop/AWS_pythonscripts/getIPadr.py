#!/usr/bin/env python
# Boto 3
import boto3
ec2 = boto3.resource('ec2')

# Boto 3
# Use the filter() method of the instances collection to retrieve
# all running EC2 instances.
instances = ec2.instances.filter(
    Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
for instance in instances:
    print(instance.id, instance.instance_type)


print(instance.public_ip_address)
print(instance.public_dns_name)

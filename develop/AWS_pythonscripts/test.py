import boto3
import boto3.ec2
help(boto3.ec2)

ec2conn = boto3.ec2.connection.EC2Connection()

reservations = ec2conn.get_all_instances()
instances = [i for r in reservations for i in r.instances]
for i in instances:
    pprint(i.__dict__)
    break # remove this to list all instances

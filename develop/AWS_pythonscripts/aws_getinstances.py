import boto3

def lambda_handler(event, context):
    ec2 = boto3.resource('ec2', region_name='eu-west-1')
    ec2_list = ec2.instances.all()

    for instance in ec2_list:
        if instance.tags:
            for tag in instance.tags:
               print(tag)
        else:
            print("no tags")


ec2 = boto3.resource('ec2', region_name='us-east-1')
ec2_list = ec2.instances.all()
print(ec2_list)
for instance in ec2_list:
    ipaddress = instance.get(u'PublicIpAddress')
    print("Ip address is ",ipaddress)
    if instance.tags:
        for tag in instance.tags:
            print(tag)
    else:
        print("no tags")



import boto3
ec2 = boto3.resource('ec2')

keyname='ec2-keypair_test2_boto3'

# create a file to store the key locally
with open(keyname+'.pem','w') as outfile:
    # call the boto ec2 function to create a key pair
    key_pair = ec2.create_key_pair(KeyName=keyname)
    # capture the key and store it in a file
    KeyPairOut = str(key_pair.key_material)
    outfile.write(KeyPairOut)


#aws_images.py
#
# Given the ID of an Amazon public AMI in one region, figure out what the
# equivalent AMI IDs are for that same AMI in all other regions known.
# If that AMI isn't defined in a region, it prints the region's name, but
# comments it out.
#

from __future__ import print_function
import boto3

# Get a starting point by finding a value in the table on this page or similar page.
# https://aws.amazon.com/amazon-linux-ami/
#
# Alternatively, you can run something like:
# aws ec2 describe-images --owners amazon --filter "Name=description,Values=Amazon Linux AMI 2017*" |
#   jq -c '.Images[] | {ImageId, Name, Description}'
#
# It will list images in the current region based on the filter expression. From there, you can
# figure out which image suits your needs. Stick the region and AMI identifier below.

# This is "amzn-ami-hvm-2017.03.1.20170623-x86_64-s3"
# This is 2017-10-03  ami-8c1be5f6 : HVM(SSD) EBS-Backed 64-bit
BASEAMIs= [ "ami-8c1be5f6"] # ,"ami-8c1be5f6", "ami-5648ad2c", "ami-cd0f5cb6" ]
BASEREGION="us-east-1"

# List all regions
client = boto3.client('ec2', region_name=BASEREGION)
regions = [region['RegionName'] for region in client.describe_regions()['Regions']]

ec2 = boto3.resource('ec2', region_name=BASEREGION)

i = 0
target_image = {}
name_filter  = {}
for aminame in BASEAMIs:
    # Figure out what the AWS Name is for this image.
    # AMI filter
    image_filter = [{'Name':'image-id', 'Values':[BASEAMIs[i]]}]
    image_names = list(ec2.images.filter(Filters=image_filter).all())
    if( len(image_names) != 1 ):
        print( "ERROR: ", len(image_names), "matched", BASEAMIs[i])
        exit( 1 )

    print( '# AMI{} is {}'.format( i, image_names[0].name ) )
    target_image[i] = image_names[0].name
    # name-based filter
    name_filter[i] = [{'Name':'name', 'Values': [target_image[i]] }]
    i = i+1

# Print the first lines of YAML
print( "Region2AMI:" )

# For every region, look up the AMI ID for that region by looking for
# the image with the same name.
for r in regions:
    ec2 = boto3.resource('ec2', region_name=r)
    i=0
    print( "  " + r + ":" )
    while (i < len(target_image.keys()) ):
        image_names = list(ec2.images.filter(Filters=name_filter[i]).all())
        if( len(image_names) != 1 ):
            print( '#  {} undefined'.format(r) )
        else:
            print( '    AMI{}: {}'.format(i, image_names[0].id ) )
        i=i+1

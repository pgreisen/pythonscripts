import boto.ec2
import sys

regions = boto.ec2.regions()


# print "Regions available: ", regions
# assert 1 == 0

for region in regions:
    try:
        ec2conn = boto.ec2.connect_to_region(region.name);
        reservations = ec2conn.get_all_reservations();

        for reservation in reservations:
                for instance in reservation.instances:
                    print instance.image_id + " " + str(instance.tags)
                    # print "Region is running: ", region

    except:
        print(sys.exc_info());



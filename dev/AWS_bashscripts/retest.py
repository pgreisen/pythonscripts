import boto.ec2
from pprint import pprint
#conn=boto.ec2.connect_to_region("eu-west-1")
conn=boto.ec2.connect_to_region("eu-central-1")
reservations = conn.get_all_instances()

print reservations

instances = [i for r in reservations for i in r.instances]

# print instances, len(instances)

# assert 1 == 0

for i in instances:
    for j in i.__dict__:
        print j
        if(j == "key_name"):

            if(i.__dict__[j] == "4dock" ):
                print i.__dict__["id"]

    # pprint(i.__dict__)

assert 1 == 0

for res in reservations:

    for inst in res.instances:

        if 'Name' in inst.tags:
            print "%s (%s) [%s]" % (inst.tags['Name'], inst.id, inst.state)
        else:
            print "%s [%s]" % (inst.id, inst.state)

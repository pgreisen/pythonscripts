# import argparse
# import googleapiclient
from googleapiclient import discovery
# import json
import googleapiclient
# from googleapiclient import discovery
# from six.moves import input
import time
# import os
# import shlex
# from google.cloud import storage
# import os
from google.cloud import *


class EnevolvGCP:

    def __init__(self):
        self.sourceimage = 'https://www.googleapis.com/compute/v1/projects/ubuntu-os-cloud/global/images/family/ubuntu-1404-lts'
        self.machine_type = "zones/%s/machineTypes/n1-standard-1"
        # docker image
        # cos-69-10895-242-0

    def create_instance(self, compute, project, zone, name, bucket):
        '''

        :param project:
        :param zone:
        :param name:
        :param bucket:
        :return:
        '''
        machine_type = self.machine_type % zone
        config = {
            'name': name,
            'machineType': machine_type,

            # Specify the boot disk and the image to use as a source.
            'disks': [
                {
                    'boot': True,
                    'autoDelete': True,
                    'initializeParams': {
                        'sourceImage': self.sourceimage
                    }
                }
            ],

            # Specify a network interface with NAT to access the public
            # internet.
            'networkInterfaces': [{
                'network': 'global/networks/default',
                'accessConfigs': [
                    {'type': 'ONE_TO_ONE_NAT', 'name': 'External NAT'}
                ]
            }]
        }
        return compute.instances().insert(
            project=project,
            zone=zone,
            body=config).execute()

    def wait_for_operation(self, compute, project, zone, operation):
        print('Waiting for operation to finish...')
        while True:
            result = compute.zoneOperations().get(
                project=project,
                zone=zone,
                operation=operation).execute()

            if result['status'] == 'DONE':
                print("done.")
                if 'error' in result:
                    raise Exception(result['error'])
                return result
            time.sleep(1)

    def list_instances(self, compute, project, zone):
        result = compute.instances().list(project=project, zone=zone).execute()
        return result['items']

    def test(self):

        compute = googleapiclient.discovery.build('compute', 'v1')
        print('Creating instance.')
        # Project
        project = 'ngsdataanalysis'
        # zone
        zone = 'us-east1-b'
        ## Name to identify instance
        name = 'ubuntu-1404-trusty-v20190429'
        ## Bucket - not defined for first run
        backet = 'empty-stuff'

        # compute.instances().start(project=project, zone=zone, instance=*).execute()
        # instance_name = 'ubuntu-1404-trusty-v20190429'
        instance_name = 'testpg'
        operation = create_instance(compute, project, zone, instance_name, bucket)

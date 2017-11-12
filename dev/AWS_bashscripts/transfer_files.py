import boto.ec2,os, shutil,subprocess
from pprint import pprint


class TransferFiles:

    def __init__(self):
        self.region = "eu-central-1"
        self.keypem = ""
        self.keypems = []
        self.ipaddress = ""
        self.path = "."
        self.pdbfile = ""
        # path to docking files in a Docking directory
        self.pathdockingfiles = "/z/insulin/users/pjug/Projects/SideProjects/GlucagonDocking/Refinement_R2_pareto/20161202_potential/files"
        self.keyname = ""

    def set_key(self):
        for file in os.listdir(self.path):
            if(file.endswith(".pem")):
                print file
                self.keypems.append(file)
            if(file.endswith(".pdb")):
                self.pdbfile = file

        assert len(self.keypems)  == 1

        self.keypem = self.keypems[0].split(".")[0]


    def setup_files(self):
        # shutil.copy(src, dst)
        copyfiles = "cp -r " + self.pathdockingfiles + "/Docking ."
        # print copyfiles
        copyrunscript = "cp -r " + self.pathdockingfiles + "/run_docking.sh ."
        # print copyrunscript
        cppdbfile = "cp "+self.pdbfile+" Docking/"
        # print cppdbfile
        # shutil.copy(path.dockingfiles, self.path)
        tardir = "tar zcf Docking.tgz Docking; rm -rf Docking;"

        ## assert 1 == 0

        subprocess.Popen(copyfiles, shell=True).wait()
        subprocess.Popen(copyrunscript, shell=True).wait()
        subprocess.Popen(cppdbfile, shell=True).wait()
        subprocess.Popen(tardir, shell=True).wait()


    def transfer_files(self):
        file2copy = "Docking.tgz "+self.pathdockingfiles+"/run_docking.sh "+self.pathdockingfiles+"/database.tgz"

        chmod_permission = "chmod 400 "+self.keypem+".pem"
        # print chmod_permission
        subprocess.Popen(chmod_permission, shell=True).wait()

        copyfiles = "scp -i "+self.keypem+".pem "+file2copy+"  ec2-user@"+self.ipaddress+":"
        # print "Copy files: ",copyfiles
        subprocess.Popen(copyfiles, shell=True).wait()



    def set_idaddress(self):
        # setup and AWS environment
        conn = boto.ec2.connect_to_region(self.region)
        reservations = conn.get_all_instances()
        # get the json format for each of the instances launched in
        # the specified region
        instances = [i for r in reservations for i in r.instances]
        # set IPaddresses
        for i in instances:
            for j in i.__dict__:
                # print j
                if(j == "key_name"):
                    if(i.__dict__[j] == self.keypem ):
                        self.ipaddress = i.__dict__["ip_address"]

    def run_script(self):

        ssh_to_instance = "ssh -i " + self.keypem + ".pem ec2-user@" + self.ipaddress

        hosts = "ssh -i " + self.keypem + ".pem ec2-user@" + self.ipaddress


        # HOST = "www.example.org"
        # Ports are handled in ~/.ssh/config since we use OpenSSH
        # COMMAND = "uname -a"
        COMMAND = "sh run_docking.sh "

        #ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
        #                       shell=False,
        #                       stdout=subprocess.PIPE,
        #                       stderr=subprocess.PIPE)

        ssh = subprocess.Popen([hosts,COMMAND],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        errdate = ssh.communicate()[1]

        # result = ssh.stdout.readlines()
        # errdata = prog.communicate()[1]


    def main(self):
        self.set_key()
        self.set_idaddress()

        print "Setting IP adress: ", self.ipaddress
        print "Setting key: ", self.keypem

        self.setup_files()
        self.transfer_files()

        # self.run_script()

        ssh_to_instance = "ssh -i " + self.keypem + ".pem ec2-user@" + self.ipaddress
        print "Debug debug ", ssh_to_instance
        subprocess.Popen(ssh_to_instance, shell=True).wait()

        # cmd = "sh run_docking.sh "
        # subprocess.Popen(cmd, shell=True).wait()



if __name__ == "__main__":
    run = TransferFiles()
    run.main()


import sys, shutil, os, subprocess, argparse

'''
Argument1 is number of dockings to setup
Argument2 is path to input pdd

'''


class SetupDocking:


    def __init__(self):
        self.exe = "/novo/users/pjug/bin/rosetta/rosetta_scripts.static.linuxgccrelease"
        self.flags = "/novo/users/pjug/templates/Docking/flags"
        self.database = "/novo/users/pjug/bin/database"
        self.xml = "/novo/users/pjug/templates/Docking/local_docking.xml"
        # full path to pdbfile
        self.pdbfile = ""
        self.number_of_jobs = 0
        # hours
        self.walltime = "10:00:00"
        # mem
        self.mem = "1gb"
        self.output_directory = ""
        self.id = ""
        self.priority = 0
        self.nstruct = "50"


    def write_to_file(self,template):
        with open("qsub.sh",'w') as  f:
            for line in template:
                f.write(line)


    def setup_directories(self):
        shutil.os.mkdir("output")
        PTH = os.path.abspath('./')
        self.output_directory = PTH+"/output/"
        for i in range(self.number_of_jobs):
            shutil.os.mkdir(str(i))
            os.chdir(str(i) )
            self.id = str(i)
            tmp = self.dump_subfile( )
            self.write_to_file(tmp)
            os.chdir("../")


    def dump_subfile(self):
        template = '''#!/bin/sh -f
#PBS -N JobName
#PBS -j oe
#PBS -l walltime='''+self.walltime+'''
#PBS -l mem='''+self.mem+'''
#PBS -m ae
#PBS -M pgreisen@gmail.com
#PBS -p '''+self.priority+'''
# priority in the queue
#  #PBS -p 0
#          Specify the number of nodes requested and the
#          number of processors per node.
#PBS -l nodes=1:ppn=1
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
'''+self.exe+''' -database '''+self.database+''' @'''+self.flags+''' -parser:protocol '''+self.xml+''' -s  '''+self.pdbfile+''' -out:nstruct '''+self.nstruct+''' -in:file:native '''+self.pdbfile+''' -out:file:silent '''+self.output_directory+str( self.id )+'''default.out -out:file:silent_struct_type binary  -out:file:scorefile '''+self.output_directory+self.id+'''dock.sc -out:prefix '''+self.id+'''>'''+self.output_directory+self.id+'''_logfile

exit'''
        return template



    # Parameters to run rosetta on hyak with patchdock output
    # locate all parameter files in one central director
    def main(self):

        parser = argparse.ArgumentParser(description="Setup jobs protein-protein docking on Davinci - location of PDB file has to be full path ")

        parser.add_argument('--xml',dest='xml', help='xml-file to be used with rosetta script',default="/novo/users/pjug/templates/Docking/local_docking.xml")

        parser.add_argument('--walltime',dest='walltime', help='Wall time limit',default='10:00:00')

        parser.add_argument('--exe',dest='exe', help='Location of rosetta script executable',default="/novo/users/pjug/bin/rosetta/rosetta_scripts.static.linuxgccrelease" )

        parser.add_argument('--flags',dest='flags', help="flags to be used for the docking",default="/novo/users/pjug/templates/Docking/flags")

        parser.add_argument('--database',dest='database', help="path to rosetta database",default="/novo/users/pjug/bin/database" )

        parser.add_argument('--pdbfile','-p',dest='pdbfile', help="Full path to pdb file" )

        parser.add_argument('-n',dest='number_of_jobs', help="Number of individual starting jobs to run ",default=10,type=int)

        parser.add_argument('--mem',dest='mem', help="Memory used for job",default="1gb" )

        parser.add_argument('--priority',dest='priority', help="Priority of jobs ( #PBS -p 0) ",default="0" )

        parser.add_argument('--nstruct',dest='nstruc', help="number of jobs per job",default="50" )


        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])

        # make sure that pdbfile exist
        ##import pdb;pdb.set_trace()
        assert os.path.isfile(self.pdbfile)

        # setup
        self.setup_directories()


        '''
        print "The following parameters have been set:"
        print "xml-file is : ", self.xml
        print "parameter file is : ", self.parameterfile
        print "setupdirs is : ", self.setupdirs
        if(self.setupdirs == "1"):
            self.setup_directories()
        PTH = os.path.abspath('./')
        outputdir = PTH+'/output/'
        try:
            shutil.os.mkdir('joblist')
            shutil.os.mkdir('output')
        except:
            print 'Joblist directory already exists'
        '''




if __name__ == "__main__":
    run = SetupDocking()
    run.main()
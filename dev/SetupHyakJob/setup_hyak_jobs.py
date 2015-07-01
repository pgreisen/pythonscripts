import sys, shutil, os, subprocess, argparse

'''
Argument1 is parameter file
Argument2 is xml-file
Argument3 is modulus 
## Not yet: Argument3 is number of processors

'''


class SetupJobsHyak:


    def __init__(self):
        self.exe = "/gscratch/baker/greisen/hyak_run/rosetta_scripts.static.linuxgccrelease"
        self.flags = "/gscratch/baker/greisen/hyak_run/flags"
        self.database = "/gscratch/baker/greisen/hyak_run/database/"
        self.setupdirs = "1"
        self.parameterfile = ""
        self.sqldatabase = "sqlname"


    def setup_directories(self):
        PTH = os.path.abspath('./')
        # Get the files in the directory
        pdbfiles = os.listdir(PTH)
        i = 1
        for pdbfile in pdbfiles:
            if( os.path.isfile(pdbfile) and str(pdbfile[-3:]) == 'pdb'):
                shutil.os.mkdir(str(i))
                shutil.move(pdbfile,str(i))

                i = i +1


    def write_subfile(self,jobname,job,outputdir,PTH,number_of_processors):
        template = '''#!/bin/bash
## RENAME FOR YOUR JOB
#PBS -N matching'''  +str(jobname)+'_'+str(job)+'''
## EDIT FOR YOUR JOB
#PBS -l walltime=72:00:00
## For 12/16 core nodes are used by the Baker lab
#PBS -l nodes=1:ppn='''+str(number_of_processors)+''',mem=22gb,feature='''+str(number_of_processors)+'''core
## EDIT FOR YOUR JOB
## Put the STDOUT and STDERR from jobs into the below directory
## 
## Put both the stderr and stdout into a single file
##PBS -j oe
##PBS -M pgreisen@gmail.com
## a mail is sent when the job is aborted by the batch system.
## b mail is sent when the job begins execution.
## e mail is sent when the job terminates.
##PBS -m abe

## Output directory
#PBS -o '''+str(outputdir)+'''

## EDIT FOR YOUR JOB
## Sepcify the working directory for this job bundle
## NOTE: we run from the LINKS directory - NOT the dir
## containing the actual script files.  This is CRITICAL
## for our checkpoinitng scheme
#PBS -d '''+str(outputdir)+'''

## If you cant run as many tasks as there are cores due to
## memory constraints you can simply set HYAK_SLOTS to a
## number instead.
#HYAK_SLOTS=6
HYAK_SLOTS=`wc -l < $PBS_NODEFILE`

## Prevent tasks from exceeding the total RAM of the node
## Requires HYAK_SLOTS to be set to number of tasks started.
NODEMEM=`grep MemTotal /proc/meminfo | awk \'{print $2}\'`
NODEFREE=$((NODEMEM-2097152))
MEMPERTASK=$((NODEFREE/HYAK_SLOTS))
ulimit -v $MEMPERTASK


module load parallel_sql
parallel_sql --sql -a parallel --exit-on-term --sql-set '''+self.sqldatabase+''' -j $HYAK_SLOTS

'''
        pbsfile = open("pbs_"+job+".sh",'w')
        pbsfile.write(template)
        pbsfile.close()
        return 1


    # Parameters to run rosetta on hyak with patchdock output
    # locate all parameter files in one central director
    def main(self):

        parser = argparse.ArgumentParser(description="Setup jobs for HYAK - for files used full path is assumed")

        parser.add_argument('--xml',dest='xml', help='xml-file to be used with rosetta script',default="")

        parser.add_argument('--np',dest='number_of_processors', help='Numbers of processors used for the simulation',default='12')

        parser.add_argument('-q',dest='queue', help='queue for the simulatio',default='bf')

        parser.add_argument('--wt',dest='wt', help='Wall time limit',default='3:59:00')

        parser.add_argument('--parameterfile',dest='parameterfile', help='Parameter file used for the design run ( full path)',default="")

        parser.add_argument('--flags',dest='flags', help='flag files',default="/gscratch/baker/greisen/hyak_run/flags")

        parser.add_argument('--database',dest='database', help='flag files',default="/gscratch/baker/greisen/hyak_run/database/")

        parser.add_argument('--setupdirs',dest='setupdirs', help='Will setup a directory full of PDB files',default="1")

        parser.add_argument('--exe',dest='exe', help='Location of rosetta script executable',default="/gscratch/baker/greisen/hyak_run/rosetta_scripts.static.linuxgccrelease")

        parser.add_argument('--sqldatabase',dest='sqldatabase', help='Name of sql database',default="sql name")

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])


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

        # Get the files in the directory
        pdbdir = os.listdir(PTH)

        with open("joblist/job",'w') as f:

            for pdbfile in pdbdir:
                if( pdbfile == 'joblist' or pdbfile == 'output'):
                    continue

                if(os.path.isdir(pdbfile) ):
                    f.write(self.exe+' -database '+self.database+' @'+self.flags+' -extra_res_fa '+self.parameterfile+' -s  '+str(PTH)+'/'+str(pdbfile)+'/*.pdb  -parser:protocol '+str(self.xml)+' > log'+str(pdbfile)+'\n')

if __name__ == "__main__":
    run = SetupJobsHyak()
    run.main()
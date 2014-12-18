import sys, shutil, os, subprocess, argparse

'''

Requires a directory with pdb files. Loops over all the pdbs and makes a directory
and submit the job.

Requires: User needs to set the path to the different parameters such as parameters, flags,
xml, database, and executables. 

and is simply executed like this:

python submit_glide.py


'''


class SubmitGlide:


    def __init__(self):
        self.xml = None
        self.params = None
        self.flags = "flags"
        self.revert_to_native = False
        self.rosetta_exe = "/work/greisen/files/glide_examples/rosetta_exe_master_03_12_2014/rosetta_scripts.static.linuxgccrelease"
        self.parameter_path = "/work/greisen/files/glide_examples/parameters/"
        self.database = "database"
        self.database_tgz = "/work/greisen/files/glide_examples/newer_database/database.tgz"
        self.directory_w_files = "/work/greisen/files/glide_examples/files/"
        self.confs = None
        self.scaffold_pdb_path = ""




    def write_wrapper(self,pdbname):
        wrapper = '''
        #!/bin/bash

        tar zxf '''+self.database+'''.tgz

        ./rosetta_scripts.static.linuxgccrelease -overwrite @'''+str(self.flags)+''' -parser:protocol '''+str(self.xml)+''' -out:prefix ${1} -extra_res_fa '''+self.params+''' -database '''+self.database+''' -s '''+str(pdbname)+''' > run_${1}.txt

        RETVAL=$?

        echo rosetta retval is $RETVAL
        if [ "$RETVAL" != "0" ]; then
        echo exiting with $RETVAL
        exit $RETVAL
        fi

        echo checking log file for "jobs considered"
        grep -q \'jobs considered\' run_${1}.txt

        RETVAL=$?
        echo grep retval is $RETVAL
    
        if [ "$RETVAL" != "0" ]; then
        echo NOT found, exiting with 1
        exit 1
        fi
    
        echo string found, exiting with 0
        exit 0
        '''
        return wrapper

    def write_wrapper_revert_to_native(self,pdbname):
        wrapper = '''
        #!/bin/bash

        tar zxf '''+self.database+'''.tgz

        ./rosetta_scripts.static.linuxgccrelease -overwrite @'''+str(self.flags)+''' -parser:protocol '''+str(self.xml)+''' -out:prefix ${1} -extra_res_fa '''+self.params+''' -database '''+self.database+''' -s '''+str(pdbname)+'''  -native native.pdb > run_${1}.txt

        RETVAL=$?

        echo rosetta retval is $RETVAL
        if [ "$RETVAL" != "0" ]; then
        echo exiting with $RETVAL
        exit $RETVAL
        fi

        echo checking log file for "jobs considered"
        grep -q \'jobs considered\' run_${1}.txt

        RETVAL=$?
        echo grep retval is $RETVAL

        if [ "$RETVAL" != "0" ]; then
        echo NOT found, exiting with 1
        exit 1
        fi

        echo string found, exiting with 0
        exit 0
        '''
        return wrapper


    def update_condor_script(self, pdbfile):
        PTH = os.path.abspath('./')
        template = '''
        notification=Never
        should_transfer_files = YES
        when_to_transfer_output = ON_EXIT
        on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)

        transfer_input_files = '''+self.database_tgz+''', '''+self.directory_w_files+self.xml+''', '''+self.directory_w_files+self.flags+''', '''+self.parameter_path+self.params+''', '''+str(PTH)+'''/'''+str(pdbfile)+''', '''+str(PTH)+'''/run_wrapper.sh, '''+self.rosetta_exe+''', '''+self.parameter_path+self.confs+'''

        Executable = run_wrapper.sh
        universe = vanilla
        copy_to_spool = false

        Error = err.$(Process)
        Output = out.$(Process)
        Log = condor_log.txt

        Arguments = $(Process)
        queue 1
    
        '''

        return template

    def update_condor_script_revert_to_native(self, pdbfile):
        PTH = os.path.abspath('./')
        template = '''
        notification=Never
        should_transfer_files = YES
        when_to_transfer_output = ON_EXIT
        on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)

        transfer_input_files = '''+self.database_tgz+''', '''+self.directory_w_files+self.xml+''', '''+self.directory_w_files+self.flags+''', '''+self.parameter_path+self.params+''', '''+str(PTH)+'''/'''+str(pdbfile)+''', '''+str(PTH)+'''/native.pdb,'''+str(PTH)+'''/run_wrapper.sh, '''+self.rosetta_exe+''', '''+self.parameter_path+self.confs+'''

        Executable = run_wrapper.sh
        universe = vanilla
        copy_to_spool = false

        Error = err.$(Process)
        Output = out.$(Process)
        Log = condor_log.txt

        Arguments = $(Process)
        queue 1

        '''

        return template



    def write_template_to_file(self,template,name='condor.submit'):
        cndr = open(name,'w')
        for line in template:
            cndr.write(line)

    def copy_native_pdb(self,pdbname):
        assert len( pdbname ) == 4
        pdbname = pdbname.lower()
        native_pdb = pdbname+'_nohet_1_relax.pdb'
        pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+native_pdb
        move_files = 'cp '+pdb_path+' .'
        subprocess.Popen(move_files,shell=True).wait()
        move_file2 = 'cp '+pdb_path+' ./native.pdb'
        subprocess.Popen(move_file2,shell=True).wait()
        return native_pdb


    def copy_native_pdb_scaffolds(self,pdbname):
        # here the pdbname is a directory
        pdbname = pdbname.upper()

        native_pdb = self.scaffold_pdb_path+pdbname+'/*_0001.pdb'

        # pdb_path = '/lab/shared/scaffolds/'+pdbname[1:3]+'/'+pdbname+'/'+native_pdb
        move_files = 'cp '+native_pdb+' .'
        subprocess.Popen(move_files,shell=True).wait()

        move_file2 = 'cp '+native_pdb+' ./native.pdb'
        subprocess.Popen(move_file2,shell=True).wait()

        return native_pdb



    def get_pdbname( self, pdbfile, split_string_number):
        #debug
        print "PDB ID is: ", pdbfile.split('_')[split_string_number]
        return  pdbfile.split('_')[split_string_number]

    def main(self):

        parser = argparse.ArgumentParser(description="Submit files for Glide - you need to specify xml-file, parameter-file ( default for flags and executable )")
        # get the initial rosetta design as input
        parser.add_argument("--xml", dest="xml", help="Name of the xml-file (located in /work/greisen/files/glide_examples/files/ ")

        parser.add_argument("-p", "--params", dest="parameterfile", help="parameterfile located in (/work/greisen/files/parameters/)")

        parser.add_argument("--revert_to_native", action="store_true", dest="revert_to_native", help="For revert to native using the native sequence this flag is necessary", default=False)

        # split pdb string to get native pdb
        parser.add_argument("-s", "--split_string_number", dest="split_string_number", help="The name of the pdb file is split by '_' and this number is used as the name for the pdb", default=3, type=int )

        parser.add_argument("-c", "--confs", dest="confs", help="The name of the ligand conformers", default=None, type=str )


        parser.add_argument("--scaffold_database", dest="scaffold_pdb_path", help="User defined scaffold set", default=None, type=str )



        input_variables = parser.parse_args()

        self.xml = input_variables.xml
        self.params = input_variables.parameterfile
        self.revert_to_native = input_variables.revert_to_native
        self.confs = input_variables.confs

        self.scaffold_pdb_path = input_variables.scaffold_pdb_path


        PTH = os.path.abspath('./')
        # Get the files in the directory
        pdbfiles = os.listdir(PTH)
        i = 1
        for pdbfile in pdbfiles:

            if(os.path.isfile(pdbfile) and str(pdbfile[-3:]) == 'pdb'):
                shutil.os.mkdir(str(i))
                shutil.move(pdbfile,str(i))
                os.chdir(str(i))
                #Condor
                condor_template = self.update_condor_script(pdbfile)
                # Wrapper
                wrapper_template = self.write_wrapper(pdbfile)

                if( self.revert_to_native ):
                    condor_template = self.update_condor_script_revert_to_native(pdbfile)
                    wrapper_template = self.write_wrapper_revert_to_native(pdbfile)


                    pdbname = self.get_pdbname( pdbfile, input_variables.split_string_number)

                    # self.copy_native_pdb( pdbname )
                    # also need a list of pdbs not present in the location.
                    if(self.scaffold_pdb_path == None):

                        self.copy_native_pdb( pdbname )


                    else:
                        self.copy_native_pdb_scaffolds( pdbname )

                        # also need a list of pdbs not present in the location.


                self.write_template_to_file(condor_template,name='condor.submit')
                self.write_template_to_file(wrapper_template,name='run_wrapper.sh')

                os.chdir('../')
                i = i +1


if __name__ == "__main__":
    run = SubmitGlide()
    run.main()

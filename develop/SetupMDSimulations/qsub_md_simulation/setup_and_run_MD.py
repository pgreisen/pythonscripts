import os,subprocess, sys
from SolvationFileForAmber import *
#from optparse import OptionParser
from CleanPDBFormat import *
from SetupSolventMinimization import *
import argparse

'''

A directory with one pdb file.

1. It will generate a holo-structure of the input pdb

2. minimization of the input structure

3. Solvation of structure

4. Minimization of solvent

5. Minimization of protein and water +ions

6. Heat up system to temperature

7. MD simulation. 


'''


class SetupAndRunMD:


    def __init__(self):
        '''


        Setting parameters necessary for running the MD simulation


        '''

        self.temperature = 300
        self.number_of_runs = 1
        self.number_of_processors = 12
        self.length_of_run = 10000000
        self.queue = "bf"
        self.wt = "3:59:00"
        self.memory = "22gb"
        self.dt = "0.002"
        self.output_freq = "5000"
        self.ff = "leaprc.ff14SB"
        self.cyclic_peptide = 0
        self.sffa = SolvationFileForAmber()



    def setup_and_run_qsub(self,wt="3:59:00",np="12",queue="bf"):
        template = '''#!/bin/bash
## RENAME FOR YOUR JOB
#PBS -N MD
## EDIT FOR YOUR JOB
#PBS -l walltime='''+wt+'''
## Cores/nodes
#PBS -l nodes=1:ppn='''+np+''',mem='''+self.memory+''',feature='''+np+'''core
## EDIT FOR YOUR JOB
## Put the STDOUT and STDERR from jobs into the below directory
## Put both the stderr and stdout into a single file
#PBS -j oe
## a mail is sent when the job is aborted by the batch system.
## b mail is sent when the job begins execution.
## e mail is sent when the job terminates.
## PBS -m abe
#PBS -d ./
mpiexec.hydra -np '''+np+''' /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i md_disulfides.in -p solv.prmtop -c heat.rst -r md.rst
'''
        tmpfile = open("qsub_md.sh",'w')
        tmpfile.write(template)
        tmpfile.close()
        exefile = open("exe.sh",'w')
        qsub_exe = "qsub -q "+queue+" qsub_md.sh"
        exefile.write(qsub_exe)
        exefile.close()
        exe = "sh exe.sh"
        subprocess.Popen(exe,shell=True).wait()

    def write_md_in(self,temperature="300",length_of_run="1500000"):
    
        template = '''
10ps MD
&cntrl
imin = 0, irest = 0, ntx = 1, iwrap=1,
ntb = 2, pres0 = 1.0, ntp = 1,
taup = 2.0,ig=-1,
cut = 10.0, ntr = 0,
ntc = 2, ntf = 2,
tempi = '''+self.temperature+''', temp0 = '''+self.temperature+''',
ntt = 3, gamma_ln = 1.0,
nstlim = '''+length_of_run+''', dt = '''+self.dt+''',
ntpr = '''+self.output_freq+''', ntwx = '''+self.output_freq+''', ntwr = '''+self.output_freq+'''
/

    '''
        md_in = open("md_disulfides.in",'w')
        md_in.write(template)
        md_in.close()

    def setup_and_run_md(self,number_of_runs,wt="3:59:00",np="12",queue="bf",temperature="300",length_of_run="5000000"):

        heat_rst_file = "cp ../HeatSystem/heat.rst ."
        heat_parm_file = "cp ../HeatSystem/*prmtop .;"

        for i in range(int( number_of_runs)):
            tmp_dir = "MD_"+str(i)
            os.mkdir( tmp_dir )
            os.chdir( tmp_dir )

            subprocess.Popen(heat_rst_file,shell=True).wait()
            subprocess.Popen(heat_parm_file,shell=True).wait()

            self.write_md_in(temperature,length_of_run)

            self.setup_and_run_qsub(wt,np,queue)
            # INSERT HERE
            os.chdir("../")

    def generate_min_sys_input_file(self):
        template = '''
Minimization of complex
 &cntrl
  imin   = 1,
  maxcyc = 20000,
  ncyc   = 10000,
  cut    = 10.0
 /
'''    
        tmp_file = open('minimize_sys.in','w')
        tmp_file.write(template)
        tmp_file.close()


    def generate_amber_exe_for_generate_parameters_for_solvation(self):
        exe = '''$AMBERHOME/bin/tleap -s -f $AMBERHOME/dat/leap/cmd/'''+self.ff+''' -f parameters_for_solvation.sh'''
        tmp_file = open('solvate_complex.sh','w')
        tmp_file.write(exe)
        tmp_file.close()

    def generate_input_file_for_library_files(self,parm,lib,ligandname):
        template='''source leaprc.gaff
loadamberparams '''+parm+'''
loadoff '''+lib+'''
OBJ = sequence { '''+ligandname+''' }
saveamberparm OBJ vacuum.prmtop vacuum.inpcrd
quit
'''
        tmp_file = open('generate_parameters.sh','w')
        tmp_file.write(template)
        tmp_file.close()


    def main(self):
        # Make the right protonation of histidines
        path = './'
        files = os.listdir(path)

        parser = argparse.ArgumentParser(description="Setup MD simulation using Amber for protein, protein-ligand or cyclic peptides")

        parser.add_argument('-t',dest='temperature', help='Temperature for the simulation',default="300")

        parser.add_argument('-n',dest='number_of_runs', help='Numbers of runs',default='1')

        parser.add_argument('--np',dest='number_of_processors', help='Numbers of processors used for the simulation',default='12')

        parser.add_argument('--lr',dest='length_of_run', help='Length of the simulation',default='1000000')

        parser.add_argument('-q',dest='queue', help='queue for the simulatio',default='bf')

        parser.add_argument('--wt',dest='wt', help='Wall time limit',default='3:59:00')

        parser.add_argument('--lib',dest='libfile', help='Name of library file',default=None)

        parser.add_argument('--parm',dest='parmfile', help='Name of parameter file',default=None)

        parser.add_argument('--lig',dest='ligname', help='Name of ligand',default=None)

        parser.add_argument('--ff',dest='ff', help='Force field used to setup the computation',default="leaprc.ff14SB")

        parser.add_argument('--cyclic_peptide',dest='cyclic_peptide', help='This will close the C-N terminal - also remember to change ff', default=0)

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])

        print "The following parameters have been set:"
        print "Temperature= ",self.temperature
        print "Number of trajectories= ",self.number_of_runs
        print "Number of processors= ",self.number_of_processors
        print "Wall time of run is= ",self.length_of_run
        print "Queue is selected= ",self.queue
        print "Parameter file is= ", self.parmfile
        print "Libfile is= ",self.libfile
        print "Ligand name is= ",self.ligname
        print "The force field is set to ", self.ff
        print "Setting AMBERHOME"

        amberhome="export AMBERHOME=/gscratch/baker/greisen/AT14/amber14"
        subprocess.Popen(amberhome,shell=True).wait()
    
        pdbfile_present = False
        pdbfile = ""
        for fl in files:
            if( os.path.isfile(fl) and fl.endswith('pdb') ):
                pdbfile = fl
                pdbfile_present = True
                print "The pdbfile is ",pdbfile
                break

        tmp = os.path.abspath("./")
        print "path is ",tmp
        # Make directory for minimization
        os.mkdir("Minimization")
        # Need to set the right libraries
        # if user has not set the right libraties
        exp1 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"
        subprocess.Popen(exp1,shell=True).wait()
        exp2="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gscratch/baker/buildbot/opt/lib/"
        subprocess.Popen(exp2,shell=True).wait()

        subprocess.Popen("mv *.pdb Minimization/",shell=True).wait()
        os.chdir("Minimization")

        ##################################################
        cpf = CleanPDBFormat()
        cpf.main(pdbfile, self.libfile, self.parmfile,self.cyclic_peptide )
        #################################################
        print "Histidines are set according to Rosettta"

        # Generate parameters for amber
        generate_parms_for_amber = "/gscratch/baker/greisen/AT14/amber14/bin/tleap -s -f /gscratch/baker/greisen/AT14/amber14/dat/leap/cmd/"+self.ff+" -f generate_input_parameters_disulfides.sh"
        subprocess.Popen(generate_parms_for_amber,shell=True).wait()

        # The script is executed in Minimization
        # Minimization of complex
        minimization = "/gscratch/baker/greisen/AT12/amber12/bin/sander -O -i /gscratch/baker/greisen/qsub_md_simulation/vacuum_min.in -p vacuum.prmtop -c vacuum.inpcrd -r min.rst"
        subprocess.Popen(minimization,shell=True).wait()

        # get minimized pdb
        get_minimized_pdb = "$AMBERHOME/bin/ambpdb -p vacuum.prmtop < min.rst > minimized.pdb"
        subprocess.Popen(get_minimized_pdb,shell=True).wait()

        os.chdir("../")

        os.mkdir("SolvateComplex")
        os.chdir("SolvateComplex")

        copy_mini_file = "cp ../Minimization/minimized.pdb ."
        subprocess.Popen(copy_mini_file,shell=True).wait()


        # copy disulfides to dir
        copy_ds_file = "cp ../Minimization/disulfide_pairs.txt ."
        subprocess.Popen(copy_ds_file,shell=True).wait()

        disulfides = self.sffa.get_disulfide_pairs()

        if( self.libfile != None  ):
            self.sffa.get_template_protein_ligand_complex(ligname,parmfile,self.libfile, disulfides)

        elif( self.cyclic_peptide != 0 ):
            self.sffa.get_template_cp( disulfides )

        else:
            self.sffa.get_template( disulfides )


        self.generate_amber_exe_for_generate_parameters_for_solvation()

        # Execute solvation script
        # Solvate the protein - here naming issues can araise
        exe_solvation = "sh solvate_complex.sh"
        subprocess.Popen(exe_solvation,shell=True).wait()

        # Initialize class
        ssm = SetupSolventMinimization()
        ssm.main("minimized.pdb")
    
        minimize_solvent_exe = "mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i minimization_of_solvent.in -p solv.prmtop -c solv.inpcrd -r solvent_minimized.rst -ref solv.inpcrd"

        subprocess.Popen(minimize_solvent_exe,shell=True).wait()

        os.chdir("../")

        os.mkdir("Minimization_whole_system")
        os.chdir("Minimization_whole_system")

        # copy file
        solv_min_rst = "cp ../SolvateComplex/solvent_minimized.rst ."
        subprocess.Popen(solv_min_rst,shell=True).wait()

        solv_parameter_file="cp ../SolvateComplex/solv.prmtop ."
        subprocess.Popen(solv_parameter_file,shell=True).wait()

        # dump parameter file for minimization of
        # whole system
        self.generate_min_sys_input_file()

        exe_minimization_whole_system="mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i minimize_sys.in -p solv.prmtop -c solvent_minimized.rst -r minimization_whole_system.rst"
        subprocess.Popen(exe_minimization_whole_system,shell=True).wait()


        os.chdir("../")

        os.mkdir("HeatSystem")
        os.chdir("HeatSystem")


        get_minimized_file ="cp ../Minimization/minimized.pdb ."
        subprocess.Popen(get_minimized_file,shell=True).wait()

        get_whole_system_rst_file="cp ../Minimization_whole_system/minimization_whole_system.rst ."
        subprocess.Popen(get_whole_system_rst_file,shell=True).wait()
        # get the parameter file
        subprocess.Popen(solv_parameter_file,shell=True).wait()

        # setup as a class instead and inherit from it
        # setup input parameter file
        heat_parameter_file_setup="/gscratch/baker/buildbot/opt/bin/python2.7 /gscratch/baker/greisen/qsub_md_simulation/setup_heat_system.py minimized.pdb "+str(self.temperature)

        subprocess.Popen(heat_parameter_file_setup,shell=True).wait()

        heat_system_exe="mpiexec.hydra -np "+self.number_of_processors+" /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i heat.in -p solv.prmtop -r heat.rst -ref minimization_whole_system.rst -c minimization_whole_system.rst"
        subprocess.Popen(heat_system_exe,shell=True).wait()

        os.chdir("../")
    
        if( self.number_of_runs > 0 ):
            self.setup_and_run_md(self.number_of_runs,self.wt,self.number_of_processors,self.queue)
    
        print "DONE"


if __name__ == "__main__":
    run = SetupAndRunMD()
    run.main()
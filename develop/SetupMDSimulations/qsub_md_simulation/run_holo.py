import os,subprocess, sys
from SolvationFileForAmber import *
from optparse import OptionParser
from CleanPDBFormat import *
from SetupSolventMinimization import *

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


class RunMD:


    def __init__(self):
        self.np = "16"

        self.amberbin = "/gscratch/baker/greisen/AT12/amber12/bin/"
        # path to executable
        self.sander_pmemd = "mpiexec.hydra -np "+self.np+"  "+self.amberbin+"pmemd.MPI "
        self.sander_mpi = "mpirun -np "+self.np+"  "+self.amberbin+"sander.MPI "
        self.sander = self.amberbin+"sander"
        self.forcefield = "$AMBERHOME/dat/leap/cmd/leaprc.ff14SB"
        self.forcefield_cp = "$AMBERHOME/dat/leap/cmd/leaprc.ff12SB_cyclic_peptide"
        self.amberhome = "export AMBERHOME=/home/ec2-user/amber14"
        self.tleap = self.amberbin+"tleap "
        self.pythonscripts = "pathtopythonscripts"

    def setup_libraries(self):
        exp1 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"
        subprocess.Popen(exp1, shell=True).wait()
        # only for Hyak?
        #exp2 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gscratch/baker/buildbot/opt/lib/"
        #subprocess.Popen(exp2, shell=True).wait()



    def setup_and_run_qsub_aws(self,np="12"):
        template = '''#
#$ -cwd
#$ -j y
#$ -pe mpi '''+str(np)+'''
#$ -S /bin/bash
#
module load openmpi-x86_64
'''+self.sander_mpi+''' -O -i md_disulfides.in -p solv.prmtop -c npt.rst -r md.rst
'''
        tmpfile = open("qsub_md.sh",'w')
        tmpfile.write(template)
        tmpfile.close()
        exefile = open("exe.sh",'w')
        qsub_exe = "qsub qsub_md.sh"
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
tempi = '''+temperature+''', temp0 = '''+temperature+''',
ntt = 3, gamma_ln = 1.0,
nstlim = '''+length_of_run+''', dt = 0.002,
ntpr = 5000, ntwx = 5000, ntwr = 5000
/

    '''
        md_in = open("md_disulfides.in",'w')
        md_in.write(template)
        md_in.close()

    def setup_and_run_md(self,number_of_runs,wt="3:59:00",np="12",queue="bf",temperature="300",length_of_run="10000000"):
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
            os.chdir("../")

    def setup_and_run_md_aws(self, number_of_runs, np="16", temperature="300",length_of_run="10000000"):
        heat_rst_file = "cp ../NPT_equilibration/npt.rst ."
        heat_parm_file = "cp ../NPT_equilibration/*prmtop .;"

        for i in range(int(number_of_runs)):
            tmp_dir = "MD_" + str(i)
            os.mkdir(tmp_dir)
            os.chdir(tmp_dir)
            subprocess.Popen(heat_rst_file, shell=True).wait()
            subprocess.Popen(heat_parm_file, shell=True).wait()
            self.write_md_in(temperature, length_of_run)
            self.setup_and_run_qsub_aws( np )
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
        exe = '''$AMBERHOME/bin/tleap -s -f '''+self.forcefield+''' -f parameters_for_solvation.sh'''
        tmp_file = open('solvate_complex.sh','w')
        tmp_file.write(exe)
        tmp_file.close()


    def generate_amber_exe_for_generate_parameters_for_solvation_cp(self):
        exe = '''$AMBERHOME/bin/tleap -s -f '''+self.forcefield_cp+''' -f parameters_for_solvation.sh'''
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

    def write_vacuum_minimization(self):
        template = '''
        Initial minimization prior to MD
 &cntrl
  imin   = 1,
  maxcyc = 10000,
  ncyc   = 5000,
  ntb    = 0,
  igb    = 0,
  cut    = 12
 /
 '''
        tmp_file = open('vacuum_min.in', 'w')
        tmp_file.write(template)
        tmp_file.close()

    def npt_equilibrate(self, temperature = "300", length_of_run = "1500000"):
        heat_rst_file = "cp ../HeatSystem/heat.rst ."
        heat_parm_file = "cp ../HeatSystem/*prmtop .;"

        tmp_dir = "NPT_equilibration"
        os.mkdir(tmp_dir)
        os.chdir(tmp_dir)
        subprocess.Popen(heat_rst_file, shell=True).wait()
        subprocess.Popen(heat_parm_file, shell=True).wait()

        template = '''
    10ps MD
    &cntrl
    imin = 0, irest = 1, ntx = 5, iwrap=1,
    ntb = 2, pres0 = 1.0, ntp = 1,
    taup = 2.0,
    cut = 10.0, ntr = 0,
    ntc = 2, ntf = 2,
    tempi = ''' + temperature + ''', temp0 = ''' + temperature + ''',
    ntt = 3, gamma_ln = 1.0,
    nstlim = 50000, dt = 0.002,
    ntpr = 5000, ntwx = 5000, ntwr = 5000
    /

        '''
        md_in = open("md_disulfides.in", 'w')
        md_in.write(template)
        md_in.close()
        npt_exe = self.sander_mpi + ''' -O -i md_disulfides.in -p solv.prmtop -c heat.rst -r npt.rst'''
        subprocess.Popen(npt_exe, shell=True).wait()
        os.chdir("../")

    def main(self):

        #
        sffa = SolvationFileForAmber()
        # Make the right protonation of histidines
        path = './'
        files = os.listdir(path)

        parser = OptionParser()
        parser.add_option('-t',dest='temperature', help='Temperature for the simulation',default="300")

        parser.add_option('-n',dest='number_of_runs', help='Numbers of runs',default='25')

        parser.add_option('--np',dest='number_of_processors', help='Numbers of processors used for the simulation',default='12')

        parser.add_option('--lr',dest='length_of_run', help='Length of the simulation',default='1000000')

        parser.add_option('-q',dest='queue', help='queue for the simulatio',default='bf')

        parser.add_option('--wt',dest='wt', help='Wall time limit',default='3:59:00')

        parser.add_option('--lib',dest='libfile', help='Name of library file',default=None)

        parser.add_option('--parm',dest='parmfile', help='Name of parameter file',default=None)

        parser.add_option('--lig',dest='ligname', help='Name of ligand',default=None)

        parser.add_option('--cyclic_peptide',dest='cyclic_peptide', help='Input is a cyclic peptide',default=False, action="store_true")


        (options,args) = parser.parse_args()

        temperature = options.temperature
        number_of_runs = options.number_of_runs
        number_of_processors = options.number_of_processors
        length_of_run = options.length_of_run
        queue = options.queue
        wt = options.wt

        libfile = options.libfile
        parmfile = options.parmfile
        ligname = options.ligname

        cyclic_peptide = options.cyclic_peptide

        print "The following parameters have been set:"
        print "Temperature= ",temperature
        print "Number of trajectories= ",number_of_runs
        print "Number of processors= ",number_of_processors
        print "Wall time of run is= ",length_of_run
        print "Queue is selected= ",queue
        print "Parameter file is= ", parmfile
        print "Libfile is= ",libfile
        print "Ligand name is= ",ligname
        print "Setting AMBERHOME"
    
        if(cyclic_peptide == True):
            print "The simulation is performed on a cyclic peptide"
        # Todo - this doesnot set variable
        amberhome = self.amberhome
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
        self.setup_libraries()

        subprocess.Popen("mv *.pdb Minimization/",shell=True).wait()
        os.chdir("Minimization")

        ##################################################
        cpf = CleanPDBFormat()
        cpf.main(pdbfile, libfile, parmfile )

        self.write_vacuum_minimization()
        #################################################
        print "Histidines are set according to Rosettta"

        # Generate parameters for amber
        generate_parms_for_amber = self.tleap+" -s -f "+self.forcefield+" -f generate_input_parameters_disulfides.sh"

        # for cyclic_peptides this has to be different
        if( cyclic_peptide == True ):
            generate_parms_for_amber = self.tleap+" -s -f "+self.forcefield_cp+" -f generate_input_parameters_disulfides.sh"


        subprocess.Popen(generate_parms_for_amber,shell=True).wait()

        # The script is executed in Minimization
        # Minimization of complex
        minimization = self.sander+" -O -i vacuum_min.in -p vacuum.prmtop -c vacuum.inpcrd -r min.rst"
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

        disulfides = sffa.get_disulfide_pairs()

        if( libfile != None  ):
            sffa.get_template_protein_ligand_complex(ligname,parmfile,libfile, disulfides)

        elif( cyclic_peptide == True ):
            sffa.get_template_cp( disulfides )

        else:
            sffa.get_template( disulfides )


        self.generate_amber_exe_for_generate_parameters_for_solvation()

        if( cyclic_peptide == True ):

            self.generate_amber_exe_for_generate_parameters_for_solvation_cp()


        # Execute solvation script
        # Solvate the protein - here naming issues can araise
        exe_solvation = "sh solvate_complex.sh"
        subprocess.Popen(exe_solvation,shell=True).wait()

        # Initialize class
        ssm = SetupSolventMinimization()
        ssm.main("minimized.pdb")
    
        minimize_solvent_exe = self.sander_mpi+" -O -i minimization_of_solvent.in -p solv.prmtop -c solv.inpcrd -r solvent_minimized.rst -ref solv.inpcrd"

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

        exe_minimization_whole_system=self.sander_mpi+" -O -i minimize_sys.in -p solv.prmtop -c solvent_minimized.rst -r minimization_whole_system.rst"
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
        # Bakerlab mainly
        # heat_parameter_file_setup="/gscratch/baker/buildbot/opt/bin/python2.7 /gscratch/baker/greisen/qsub_md_simulation/setup_heat_system.py minimized.pdb "+str(temperature)

        heat_parameter_file_setup = "python "+self.pythonscripts+"setup_heat_system.py minimized.pdb " + str(temperature)
        subprocess.Popen(heat_parameter_file_setup,shell=True).wait()

        heat_system_exe=self.sander_mpi+" -O -i heat.in -p solv.prmtop -r heat.rst -ref minimization_whole_system.rst -c minimization_whole_system.rst"
        subprocess.Popen(heat_system_exe,shell=True).wait()

        os.chdir("../")

        # Equilibrate 100 ps in the NPT ensemble


        self.npt_equilibrate()

        if( number_of_runs > 0 ):
            self.setup_and_run_md_aws(number_of_runs,number_of_processors)
        print "DONE"


if __name__ == "__main__":
   run = RunMD()
   run.main()
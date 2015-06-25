import os,subprocess, sys
from SolvationFileForAmber import *
from optparse import OptionParser



'''

A directory with one pdb file.

1. It will generate a apo-structure of the input pdb

2. minimization of the input structure

3. Solvation of structure

4. Minimization of solvent

5. Minimization of protein and water +ions

6. Heat up system to temperature

7. MD simulation. 


'''
sffa = SolvationFileForAmber()

def setup_and_run_qsub(wt="3:59:00",np="12",queue="bf"):
    template = '''#!/bin/bash
## RENAME FOR YOUR JOB
#PBS -N MD
## EDIT FOR YOUR JOB
#PBS -l walltime='''+wt+'''
## For 12 core nodes are used by the Baker lab
#PBS -l nodes=1:ppn='''+np+''',mem=22gb,feature='''+np+'''core
## EDIT FOR YOUR JOB
## Put the STDOUT and STDERR from jobs into the below directory
## 
## Put both the stderr and stdout into a single file
#PBS -j oe
##PBS -M pgreisen@gmail.com
## a mail is sent when the job is aborted by the batch system.
## b mail is sent when the job begins execution.
## e mail is sent when the job terminates.
## PBS -m abe

## Output directory
## Output directory

## EDIT FOR YOUR JOB
## Sepcify the working directory for this job bundle
## NOTE: we run from the LINKS directory - NOT the dir
## containing the actual script files.  This is CRITICAL
## for our checkpoinitng scheme
#PBS -d ./

mpiexec.hydra -np '''+np+''' /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i md_disulfides.in -p apo_solv.prmtop -c heat_apo.rst -r md.rst
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

def write_md_in(temperature="300",length_of_run="1500000"):
    
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

def setup_and_run_md(number_of_runs,wt="3:59:00",np="12",queue="bf"):

    heat_rst_file = "cp ../HeatSystem/heat_apo.rst ."
    heat_parm_file = "cp ../HeatSystem/*prmtop .;"

    for i in range(int( number_of_runs)):
        tmp_dir = "MD_"+str(i)
        os.mkdir( tmp_dir )
        os.chdir( tmp_dir )

        subprocess.Popen(heat_rst_file,shell=True).wait()
        subprocess.Popen(heat_parm_file,shell=True).wait()

        write_md_in(temperature,length_of_run)

        setup_and_run_qsub(wt,np,queue)
        # INSERT HERE
        os.chdir("../")

def generate_min_sys_input_file():
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
    #for line in template:
    tmp_file.write(template)
    tmp_file.close()



def generate_amber_exe_for_generate_parameters_for_solvation():
    exe = '''$AMBERHOME/bin/tleap -s -f $AMBERHOME/dat/leap/cmd/leaprc.ff12SB -f parameters_for_solvation_apo.sh'''
    tmp_file = open('solvate_complex_apo.sh','w')
    tmp_file.write(exe)
    tmp_file.close()



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

(options,args) = parser.parse_args()

temperature = options.temperature
number_of_runs = options.number_of_runs
number_of_processors = options.number_of_processors
length_of_run = options.length_of_run
queue = options.queue
wt = options.wt

print temperature
print number_of_runs
print number_of_processors
print length_of_run
print queue

print "Setting AMBERHOME"

amberhome="export AMBERHOME=/gscratch/baker/greisen/AT12/amber12"
subprocess.Popen(amberhome,shell=True).wait()

pdbfile_present = False

for fl in files:
    if( os.path.isfile(fl) and fl.endswith('pdb') ):
        pdbfile = fl
        pdbfile_present = True
        break

assert pdbfile_present == True
print "The pdbfile is ",pdbfile
tmp = os.path.abspath("./")
print "path is ",tmp

# Need to set the right libraries
# if user has not set the right libraties
exp1 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"
subprocess.Popen(exp1,shell=True).wait()
exp2="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gscratch/baker/buildbot/opt/lib/"
subprocess.Popen(exp2,shell=True).wait()

keep_histidine_as_rosetta_made_them = '/gscratch/baker/buildbot/opt/bin/python2.7 /gscratch/baker/greisen/qsub_md_simulation/change_ligand_names_and_remove_protein_hydrogens_and_insert_disulfides.py '+str(pdbfile)
subprocess.Popen(keep_histidine_as_rosetta_made_them,shell=True).wait()

print "Histidines are set according to Rosettta"
#######################################################################
# Create directory Minimization and move file to this directory
os.mkdir("Minimization")
subprocess.Popen("mv *.pdb Minimization/",shell=True).wait()
os.chdir("Minimization")

# Make sure only protein atoms are present
get_protein_atoms = "grep '^ATOM*' *.pdb  >> design1.pdb"
subprocess.Popen(get_protein_atoms,shell=True).wait()

# This file generate_input_parameters_disulfides.sh is generate when the 
# names are changed and it dumps a file which is used by amber tools
# The file was dump in the main directory hence it is moved
subprocess.Popen("mv ../generate_input_parameters_disulfides.sh ./",shell=True).wait()


# Generate parameters for amber
generate_parms_for_amber = "/gscratch/baker/greisen/AT12/amber12/bin/tleap -s -f /gscratch/baker/greisen/AT12/amber12/dat/leap/cmd/leaprc.ff12SB -f generate_input_parameters_disulfides.sh"
subprocess.Popen(generate_parms_for_amber,shell=True).wait()


# The script is executed in Minimization
# Minimization of complex
minimization = "/gscratch/baker/greisen/AT12/amber12/bin/sander -O -i /gscratch/baker/greisen/qsub_md_simulation/vacuum_min.in -p apo_vacuum.prmtop -c apo_vacuum.inpcrd -r min.rst"
subprocess.Popen(minimization,shell=True).wait()

# get minimized pdb
get_minimized_pdb = "$AMBERHOME/bin/ambpdb -p apo_vacuum.prmtop < min.rst > minimized.pdb"
subprocess.Popen(get_minimized_pdb,shell=True).wait()

os.chdir("../")
assert 1 == 0
os.mkdir("SolvateComplex")
os.chdir("SolvateComplex")

copy_mini_file = "cp ../Minimization/minimized.pdb ."
subprocess.Popen(copy_mini_file,shell=True).wait()

# copy disulfides to dir
copy_ds_file = "cp ../disulfide_pairs.txt ."
subprocess.Popen(copy_ds_file,shell=True).wait()

generate_amber_exe_for_generate_parameters_for_solvation()

sffa.main()

# Execute solvation script
# Solvate the protein - here naming issues can araise
exe_solvation = "sh solvate_complex_apo.sh"
subprocess.Popen(exe_solvation,shell=True).wait()

# Copy parameter file
copy_parameter_file="/gscratch/baker/buildbot/opt/bin/python2.7 /gscratch/baker/greisen/qsub_md_simulation/setup_solvent_minimization.py minimized.pdb"
subprocess.Popen(copy_parameter_file,shell=True).wait()

minimize_solvent_exe = "mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i minimization_of_solvent.in -p apo_solv.prmtop -c apo_solv.inpcrd -r solvent_minimized_apo.rst -ref apo_solv.inpcrd"

subprocess.Popen(minimize_solvent_exe,shell=True).wait()

os.chdir("../")

os.mkdir("Minimization_whole_system")
os.chdir("Minimization_whole_system")

# copy file
solv_min_rst = "cp ../SolvateComplex/solvent_minimized_apo.rst ."
subprocess.Popen(solv_min_rst,shell=True).wait()

solv_parameter_file="cp ../SolvateComplex/apo_solv.prmtop ."
subprocess.Popen(solv_parameter_file,shell=True).wait()

# dump parameter file for minimization of
# whole system
generate_min_sys_input_file()

exe_minimization_whole_system="mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i minimize_sys.in -p apo_solv.prmtop -c solvent_minimized_apo.rst -r minimization_whole_system_apo.rst"
subprocess.Popen(exe_minimization_whole_system,shell=True).wait()


os.chdir("../")

os.mkdir("HeatSystem")
os.chdir("HeatSystem")


get_minimized_file ="cp ../Minimization/minimized.pdb ."
subprocess.Popen(get_minimized_file,shell=True).wait()

get_whole_system_rst_file="cp ../Minimization_whole_system/minimization_whole_system_apo.rst ."
subprocess.Popen(get_whole_system_rst_file,shell=True).wait()
# get the parameter file
subprocess.Popen(solv_parameter_file,shell=True).wait()

# setup as a class instead and inherit from it
# setup input parameter file
heat_parameter_file_setup="/gscratch/baker/buildbot/opt/bin/python2.7 /gscratch/baker/greisen/qsub_md_simulation/setup_heat_system.py minimized.pdb "+str(temperature)

subprocess.Popen(heat_parameter_file_setup,shell=True).wait()

heat_system_exe="mpiexec.hydra -np 12 /gscratch/baker/greisen/AT12/amber12/bin/pmemd.MPI -O -i heat.in -p apo_solv.prmtop -r heat_apo.rst -ref minimization_whole_system_apo.rst -c minimization_whole_system_apo.rst"
subprocess.Popen(heat_system_exe,shell=True).wait()

os.chdir("../")

setup_and_run_md(number_of_runs,wt,number_of_processors,queue)

print "DONE"

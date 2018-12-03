from subprocess import call
import glob
import subprocess
dst_result = "results.txt"
subprocess.Popen("rm "+dst_result, shell=True, stdout=subprocess.PIPE).wait()
for i in glob.glob("*.pdb"):
    exe = "python ~/pythonscripts/develop/analyse_designs/AnalyseMutations.py -n ../1.pdb -d "+i+" >> "+dst_result;

    subprocess.Popen(exe, shell=True, stdout=subprocess.PIPE).wait()    
    pymol_name = i.split('.')[0]
    subprocess.Popen("mv mut.pml "+pymol_name+".pml" , shell=True, stdout=subprocess.PIPE).wait()

unique_mutations = {}

with open(dst_result,'r') as f:
    for line in f:
        pdbfile,mutations = line.split(':')
        unique_mutations[mutations] = pdbfile
dst_pml="analyse_designs.pml"

modulus=10
dummy=1
prefix=str(dummy)
for j in unique_mutations.keys():
    if(dummy % modulus == 0):
        prefix = str(dummy)
    cmd_pml = "cat "+unique_mutations[j].replace('.pdb','.pml')+" >> "+prefix+"_"+dst_pml+"; echo >> "+prefix+"_"+dst_pml;
    subprocess.Popen(cmd_pml, shell=True, stdout=subprocess.PIPE).wait()
    dummy += 1
        
    
    

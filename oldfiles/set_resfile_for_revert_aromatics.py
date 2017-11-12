#!/usr/bin/env python
import sys,subprocess
import os,re
from numpy import *
'''
The algorithm is pretty dependent on the filenames at the moment
so it will probably crash used by others at the moment.

should be run after revert to native.

'''

def get_file(pdbfile):
    pdb = []
    fl = open(pdbfile,'r')
    for line in fl:
        pdb.append(line)
    fl.close()
    return pdb


# Requires a pdb file
# Returns sequence 
def get_seq(filename):
    tmp = open(filename,'r')
    resid = ''
    seq = []
    j = 0
    for line in tmp:
        tm = line.split()
        if tm[0] == 'ATOM' and j == 0:
            resid = tm[5]
            seq.append(tm[3]+' '+tm[5])
        j = j+1
        if tm[0] == 'ATOM' and tm[5] != resid:
            seq.append(tm[3]+' '+tm[5])
            resid = tm[5]
    return seq


def get_seq_pdbfile(pdbfile):
    resid = '0'
    seq = []
    j = 0
    for line in pdbfile:
        tm = line.split()
        if line[0:4] == 'ATOM' and resid != line[22:26]: #resid:
            resid = line[22:26]
            seq.append(line[17:20]+' '+tm[5])
    return seq



# Requires two arrays (same length right now)
# Wild type = a and b = rosetta output
# Returns the a list of the mutant's mutations
def get_def(a,b):
    aromatics = ['TRP','PHE','TYR']
    # Mutated positions
    mt = []
    # Native aromatic position which are not aromatic in design
    n_aro = []
    ln = len(a)
    i = 0
    mut = 0
    for i in range(ln):
        if a[i] != b[i]:
            mut = mut +1
            mt.append(b[i].split()[1])
            if a[i].split()[0] in aromatics and b[i].split()[0] not in aromatics:
                n_aro.append(a[i].split()[1])
            
    return mut,n_aro,mt



def get_residue_energy(pdbfile,residue_energy,energy_term):

    for line in pdbfile:
        tmp_eval = line.replace('_',' ').split()

        if line[0:3] in residue_energy:

            
            if(line[3:11] == '_p:Nterm' or line[3:11] == '_p:Cterm'):
                #print line
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])


            elif(line[0:9] == 'HIS_D_p:N' or line[0:9] == 'HIS_D_p:C'):
                    #print tmp_eval
                    residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+3])
                
            elif(line[3:5] == '_p'):
                
                # print tmp_eval
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])


            # A histidine fix
            elif(line[0:5] == 'HIS_D'):
                # print tmp_eval
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])

            else:
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+1])


    return residue_energy



def get_residue_energy(pdbfile):

    residue_energy = {'TRP' : -6.5,
                      'PHE' : -5.33,
                      'TYR' : -5.68}

    positions = {'TRP' : [],
                 'PHE' : [],
                 'TYR' : []}
    
    for line in pdbfile:
        
        tmp_eval = line.replace('_',' ').split()

        if tmp_eval[0] in residue_energy:

            if(line[3:5] == '_p'):
                print tmp_eval
                if float(tmp_eval[3]) > residue_energy[tmp_eval[0]]:
                    positions[tmp_eval[0]].append(tmp_eval[2])

            else:
                if float(tmp_eval[2]) > residue_energy[tmp_eval[0]]:
                    positions[tmp_eval[0]].append(tmp_eval[1])

    
    return positions


# @require pdbfile and atom specifier
# @return list with protein or ligand

def get_ligand_coordinates(pdbfile,chemical='HETA'):
    ligand_coor = []
    for line in pdbfile:
        if line[0:4] == chemical:
            x = str(line[30:38]).rstrip()
            y = str(line[38:46]).rstrip()
            z = str(line[46:54]).rstrip()
            ligand_coor.append(array([float(x),float(y),float(z)]))
    return ligand_coor


def get_length_protein(pdbfile):
    first = '0'
    start = ''
    end = ''
    for line in pdbfile:
        if line[0:4] == 'ATOM':
            if first == '0':
                start = str(line[23:26]).rstrip()
                first = '1'
            elif first == '1':
                end = str(line[23:26]).rstrip()
    return int(start),int(end)




def main():
    DISTANCE1 = 6
    DISTANCE2 = 8
    
    path = './'
    files = os.listdir(path)
    output_string = ''
    for fl in files:
        evl = fl.find('_in.pdb')
        # evl2 = re.match('(.*)_in_00([0-9]*).pdb',fl)
        evl3 = re.match('(.*)_in_00([0-9]*)_00([0-9]*).pdb',fl)
        if(evl != -1):
            native = fl
        #    print 'The native pdb is: ',native
        # if(evl2):
        #    design = fl
        #    print 'The design file is: ',design
        if(evl3):
            revert = fl
            new_file_name = 'cp '+fl+' revert_aromatic.pdb'
            subprocess.Popen(new_file_name,shell=True).wait()
        #    print 'The revert file is: ',revert

    native = get_file(native)

    a = get_seq_pdbfile(native)
    # b = get_seq_pdbfile(design)
    revert = get_file(revert)
    c = get_seq_pdbfile(revert)

    # native_design = get_def(a,b)
    native_revert,native_aromatic_positions,positions_mutated = get_def(a,c)
    # print 'Native position',positions_mutated

    # Getting the residue energy from the reverted structure
    positions_with_bad_energy_scores = get_residue_energy(revert)



    # Next generate resfile
    # 1. get ligand coordinates
    ligand = get_ligand_coordinates(revert)


    # Collecting residue id for residues within or outside
    # ligand
    nataa = []
    natro = []
    design = []

    # Get list of pdb-file
    # pdbfl = get_pdbfile(pdbfile)

    
    for i in ligand:

        for line in revert:

            if line[0:4] == 'ATOM':
                
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()

                tmp_vector = array([float(x),float(y),float(z)])

                tmp_length = linalg.norm(tmp_vector - i)

                pdb_resid = int(str(line[23:26]).rstrip())
                    
                
                if tmp_length > DISTANCE1 and tmp_length < DISTANCE2:
                    
                    pdb_resid = int(str(line[23:26]).rstrip())
                    nataa.append(pdb_resid)


                elif tmp_length > DISTANCE2:

                    pdb_resid = int(str(line[23:26]).rstrip())
                    natro.append(pdb_resid)
                    
                    
                elif tmp_length < DISTANCE1:

                    pdb_resid = int(str(line[23:26]).rstrip())
                    design.append(pdb_resid)

                else:
                    print 'Bug in program'




    trp = []
    phe = []
    tyr = []

    
    for i in positions_with_bad_energy_scores:
        if(i == 'TRP'):
            #           if(len(positions_with_bad_energy_scores[i]) > 0):
            for j in positions_with_bad_energy_scores[i]:
                if j in positions_mutated:
                    trp.append(float(j))
        if(i == 'TYR'):
            #            if(len(positions_with_bad_energy_scores[i]) > 0):
            for k in positions_with_bad_energy_scores[i]:
                if k in positions_mutated:
                    tyr.append(float(k))
                    
        if(i == 'PHE'):
            for l in positions_with_bad_energy_scores[i]:
                if l in positions_mutated:
                    phe.append(float(l))




    # get the rigth amino acids
    nataa = set(nataa).difference(set(trp))
    nataa = set(nataa).difference(set(phe))
    nataa = set(nataa).difference(set(tyr))

    natro = set(natro).difference(set(nataa))

    natro = set(natro).difference(set(trp))
    natro = set(natro).difference(set(tyr))
    natro = set(natro).difference(set(phe))

    



    start,end = get_length_protein(revert)

    # resfile
    rs_file = open('resfile','w')
    rs_file.write('start\n')
    
    rs_1 = ' A NATAA\n'
    rs_2 = ' A NATRO\n'


    rs_trp = ' A NOTAA W\n'
    rs_tyr = ' A NOTAA Y\n'
    rs_phe = ' A NOTAA F\n'


    
    while start <= end:
        
        if start in natro:
            rs_file.write('\t'+str(start)+rs_2)

        elif start in nataa:
            rs_file.write('\t'+str(start)+rs_1)

        elif start in trp:
            rs_file.write('\t'+str(start)+rs_trp)

        elif start in tyr:
            rs_file.write('\t'+str(start)+rs_tyr)

        elif start in phe:
            rs_file.write('\t'+str(start)+rs_phe)

        start = start + 1


if __name__ == "__main__":
    main()

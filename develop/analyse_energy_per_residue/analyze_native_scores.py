#!/usr/bin/env python2.7
import sys,subprocess
import os,re
from numpy import *
from pylab import hist,savefig,clf,xlim

'''


'''

def get_file(pdbfile):
    '''
    Reads the pdb file

    '''

    pdb = []
    fl = open(pdbfile,'r')
    for line in fl:
        pdb.append(line)
    fl.close()
    return pdb



def get_residue_energy(pdbfile,residue_energy,energy_term,fl):
    '''


    '''

    for line in pdbfile:

        tmp_eval = line.replace('_',' ').split()

        if line[0:3] in residue_energy:

            
            if(line[3:11] == '_p:Nterm' or line[3:11] == '_p:Cterm'):

                if(float(tmp_eval[energy_term+2]) > 20):

                    print 'The energy term is quite high',tmp_eval[energy_term+2],tmp_eval[0],fl 

                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])


            elif(line[0:9] == 'HIS_D_p:N' or line[0:9] == 'HIS_D_p:C'):

                if(float(tmp_eval[energy_term+3]) > 20):

                    print 'The energy term is quite high',tmp_eval[energy_term+3],tmp_eval[0],fl  


                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+3])
                
            elif(line[3:5] == '_p'):


                if(float(tmp_eval[energy_term+2]) > 20):
                    print 'The energy term is quite high',tmp_eval[energy_term+2],tmp_eval[0],fl 

                
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])


            # A histidine fix
            elif(line[0:5] == 'HIS_D'):

                if(float(tmp_eval[energy_term+2]) > 20):
                    print 'The energy term is quite high',tmp_eval[energy_term+2],tmp_eval[0],fl 


                
                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+2])

            else:

                if(float(tmp_eval[energy_term+1]) > 20):
                    print 'The energy term is quite high',tmp_eval[energy_term+1],tmp_eval[0],fl 

                residue_energy[tmp_eval[0]].append(tmp_eval[energy_term+1])

    return residue_energy


def get_position_energy_term(pdbfile,energy_term):
    '''
    Get the position of the energy term e.g., total_energy which part of the line should 
    it get. 
    '''
    energy_position = 0
    for line in pdbfile:
        tmp = line.split()

        if tmp[0] == 'label':
            for j in tmp:
                if(j == energy_term):
                    break
                energy_position = energy_position + 1
    return energy_position

def write_to_file(residue_energy,energyterm):
    tmp_string = str(energyterm)+'.dat'
    tmp_file = open(tmp_string,'w')
    for key in residue_energy:
        tmp_file.write(key+' : ')
        for i in residue_energy[key]:
            tmp_file.write(i+' ')
        tmp_file.write('\n')
    tmp_file.close()

def write_data_per_residue(name,data):
    dat_per_residue = open(name+'.dat','w')
    for i in data:
        dat_per_residue.write(str(i)+'\n')


def plot_results_per_residue(residue_energy):

    statistics_file = open('stat_summary.txt','w')
    
    statistics_file.write('Residue\tMean\tS.D\tMin\tMax\tNumber\n')
    
    for key in residue_energy:
        
        tmp_energy = [ float(y) for y in residue_energy[key] ]

        write_data_per_residue(key,tmp_energy)
        
        hist(tmp_energy,bins=50)
        xlim([min(tmp_energy), max(tmp_energy)])
        savefig(key+'.png')
        clf()

        statistics_file.write(key+'\t'+str(round(mean(tmp_energy),2))+'\t'+str(round(sqrt(var(tmp_energy)),2))+'\t'+str(round(min(tmp_energy),2))+'\t'+str(round(max(tmp_energy),2))+'\t'+str(len(tmp_energy))+'\n')


def main():

    residue_energy = {
        'ALA' : [],
        'ARG' : [],
        'ASN' : [],
        'ASP' : [],
        'CYS' : [],
        'GLU' : [],
        'GLN' : [],
        'GLY' : [],
        'HIS' : [],
        'ILE' : [],
        'LEU' : [],
        'LYS' : [],
        'MET' : [],
        'PRO' : [],
        'SER' : [],
        'PHE' : [],
        'THR' : [],
        'TRP' : [],
        'TYR' : [],
        'VAL' : []
        }
       
    # term_to_analyse = ['total','fa_atr','fa_sol','fa_dun','fa_rep']
    # term_to_analyse = ['fa_rep']
    # term_to_analyse = ['fa_atr','fa_sol','fa_dun','fa_rep']
    term_to_analyse = ['total']
    path = './'

    
    # We loop over all the directories in the main 
    # directory
    dirs = os.listdir(path)

    # The energy term to be analysed
    for j in term_to_analyse:
        
        for dr in dirs:
            if os.path.isdir(dr):
                # We change directory 
                os.chdir(dr)
                
                files = os.listdir(path)
                # Loop over all the file in the directory
                for fl in files:
                    # This file should contain the scores
                    evl = fl.find('relax_0001.pdb')

                    if (evl != -1):
                        
                        pdbfile = get_file(fl)
                        
                        position_of_energy_term = get_position_energy_term(pdbfile,j)

                        residue_energy = get_residue_energy(pdbfile,residue_energy,position_of_energy_term,fl)

                os.chdir('../')

        write_to_file(residue_energy,j)
    plot_results_per_residue(residue_energy)


if __name__ == "__main__":
    main()

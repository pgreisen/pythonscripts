#!/usr/bin/env python
from os import system,popen
import string
from sys import argv
from translate_rotate import *
import sys
import math
from math import sin,cos
from numpy import *
import numpy as np
from numpy import linalg as LA

import argparse
from xyzMath import *


class ParametricHelicalBundle:


    def __init__(self):
        self.deg_to_rad= math.pi/180.
        self.R1 = 2.26
        self.w1 = 720./7.
        self.atom_num=1
        self.res_num=1

        # parameter from the parameter file:
        self.Nres = 0
        self.num_chain = 0
        self.num_to_output = 0
        self.orientation = 0
        self.R0 = 0
        self.w0 = 0
        self.ph = 0
        self.z = 0
        self.output_file_prefix = 0
        self.chain_name = 0
        self.chain_order = 0

        self.chain_set=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

        self.debug = True


    def get_translation_and_rotation(self,a=None,b=None,c=None):

        translation = a
        e1 = (a-b) / LA.norm(a-b)

        tmp_e3 = np.cross(e1,c-b)
        e3 = tmp_e3 / LA.norm( tmp_e3 )

        tmp_e2 = np.cross(e1,e3)
        e2 = tmp_e2 / LA.norm( tmp_e2 )
        print e1, e2, e3

        assert LA.norm(e1) == LA.norm(e2) == LA.norm(e3) == 1

        rotation_matrix = ([e1[0],e1[1],e1[2]],[e2[0],e2[1],e2[2]],[e3[0],e3[1],e3[2]])

        return translation, rotation_matrix



    def enumerate_combinations(self, list):

        num = len(list)
        ncomb=1
        counter=[]
        for i in range(num):
            ncomb=ncomb*list[i]
            counter.append(0)
        add = 1

        combs=[]
        for j in range(ncomb):
            combs.append(string.join(map(lambda x: str(counter[x]), range(num))))

        for i in range(num):
            counter[i]=counter[i]+add
            if counter[i]==list[i]:
                counter[i]=0
            else:
                break
        return(combs)

    def Generate_Pdbs(self, Nres,num_chain,num_to_output,w0,R,orientation,helix_phase,delta_z,output_file_name,chain_name, chain_order):

        # chain parameters
        ph = 360/num_chain
        phase=[]
        for i in range(num_chain):
            phase.append(i*ph)

        # assign a chain id to each chain
        chain_num=self.chain_set[0:num_chain]

        #rise per residue d  fixed. this constrains pitch (alpha)
        d=1.51
        z1=0.
        z2=0.
        line='ATOM      7  CA  GLY A   2       5.520   2.352   1.361  1.00 20.00     '
        last=line[54:-1]

        Res_id=[]
        CA_list=[]
        alpha=math.asin(R*w0*self.deg_to_rad/d)
        for iter in range(num_to_output):
            CA_list.append([])
            Res_id.append([])
            chain = chain_name[iter]
            orient = orientation[iter]

            if orient == 1:
                res_num=0
            else:
                res_num=Nres+1

            supercoil_phase=phase[iter]+delta_z[iter]*math.tan(alpha)/(R*self.deg_to_rad)
            # Loop over all the residue specified in the input file
            for t in range(Nres+2):      ## need two extra residues to guide placement of 1st and last residue

                a0=(w0*t+supercoil_phase)*self.deg_to_rad
                # PG note if we need to subsample the value of w1+w0 = CONSTANT Here we need a
                # modification

                # w1, t
                a1=(self.w1*t+helix_phase[iter]+phase[iter])*self.deg_to_rad   # set ref point for phase to be along supercoil radius

                x=R*math.cos(a0) + self.R1*math.cos(a0) * cos(a1) - self.R1*cos(alpha)*sin(a0)*sin(a1)
                y=R*math.sin(a0) + self.R1*sin(a0)*cos(a1) + self.R1*cos(alpha)*cos(a0)*sin(a1)
                z= R*w0*t*self.deg_to_rad/math.tan(alpha +.00000001)-self.R1*sin(alpha)*sin(a1)+delta_z[iter]

                CA_list[iter].append( (res_num,array([x,y,z])) )

                Res_id[iter].append(res_num)
                self.atom_num=self.atom_num+1
    
                if orient  == 1:
                    res_num=res_num+1
                else:
                    res_num=res_num-1

        # convert CA trace to full backbone model by superimposing on ideal template
        # by matching 3 consecutive CA atoms
        # set up ideal template
        stub_file=map(string.split,open('ideal.pdb','r').readlines())
        atom=[]

        tmp_atom = []
        for line in stub_file:
            atom.append( (array([float(line[6]),float(line[7]),float(line[8] ) ])))
            tmp_atom.append( Vec( float(line[6]), float(line[7]) , float(line[8]) ) )

        # rotation and translational vectors
        t_v, r_m = self.get_translation_and_rotation( atom[6], atom[1], atom[11]  )

        if(self.debug == True):

            print "Vec1",tmp_atom[6]
            print "Vec2",tmp_atom[1]
            print "Vec3",tmp_atom[11]

        ideal_stub=stub( tmp_atom[6], tmp_atom[1], tmp_atom[11] )

        if(self.debug == True):

            print "The ideal stub is: ", ideal_stub

        full_pdb=open(output_file_name,'w')
        # this value is now set in the constructor
        atom_num=1

        print "Length of CA_list", len(CA_list)

        res_num=0
        for counter in range(num_to_output):
            iter=int(chain_order[counter])
            chain=chain_name[iter]
            CA_chain_u=CA_list[iter]
            CA_chain = sorted(CA_chain_u, key = lambda res: res[0])
            for res in range(1,Nres+1):
                res_num=res_num+1

                # STUB INSERT ROTATION/TRANSLATION
                ##actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])

                tmp_vec_a = Vec( float(CA_chain[res][1][0]),float(CA_chain[res][1][1]),float(CA_chain[res][1][2]) )
                tmp_vec_b = Vec( float(CA_chain[res-1][1][0]),float(CA_chain[res-1][1][1]),float(CA_chain[res-1][1][2]) )
                tmp_vec_c = Vec( float(CA_chain[res+1][1][0]),float(CA_chain[res+1][1][1]),float(CA_chain[res+1][1][2]) )

                print "tmp_vec_a", tmp_vec_a
                print "tmp_vec_b", tmp_vec_b
                print "tmp_vec_c", tmp_vec_c


                old_actual_stub = stub(tmp_vec_a,tmp_vec_b, tmp_vec_c)

                ##actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])


                #transform=actual_stub * ~ideal_stub
                transform=old_actual_stub * ~ideal_stub

                print "The transform is equal to ", transform

                coords=transform*atom[5]
                full_pdb.write('ATOM %6d  N   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
                self.atom_num=self.atom_num+1

                # CA   (use actual CA from trace rather than superimposed one)
                coords=CA_chain[res][1]
                tcoords=transform*atom[6]

                full_pdb.write('ATOM %6d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
                self.atom_num=self.atom_num+1

                #  NH
                coords=transform*atom[7]
                full_pdb.write('ATOM %6d  H   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
                self.atom_num=self.atom_num+1

                #  C
                coords=transform*atom[8]
                full_pdb.write('ATOM %6d  C   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
                self.atom_num=self.atom_num+1

                # O
                coords=transform*atom[9]
                full_pdb.write('ATOM %6d  O   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
                self.atom_num=self.atom_num+1

                start_d=Vec.distance(atom[8],atom[6])
                end_d = Vec.distance(transform*atom[8],transform*atom[6])


    def input_params(self,input):

        output_file_prefix=input[0][0]
        Nres = int(input[1][0])
        w0_begin = float(input[2][0])
        w0_end = float(input[2][1])
        w0_iter=int(input[2][2])
        R0_begin  = float(input[3][0])
        R0_end = float(input[3][1])
        R0_iter = int(input[3][2])
        num_chain=int(input[4][0])
        num_to_output=int(input[5][0])
        orientation=[]
        for i in range(num_to_output):
            orientation.append( int(input[6][i] ))
        phase=[]
        for i in range(num_to_output):   ## restrict phase and delta_z to switch sign in sym_mates
            phase.append( (float(input[7+i][0]),float(input[7+i][1]),int(input[7+i][2])))
        z_list=[]
        for i in range(num_to_output):
            l=7+num_to_output+i
            z_list.append( ( float(input[l][0]),float(input[l][1]),int(input[l][2]) ))
        chain_name=[]
        for i in range(num_to_output):
            chain_name.append(input[l+1][i])
        chain_order=[]
        for i in range(num_to_output):
            chain_order.append(input[l+2][i])
        
        print ' ###############    SUPERHELIX PARAMS   ############## \n'
        print ' output file prefix:  %s \n'%(output_file_prefix)
        print ' helix length:  %s \n'%(Nres)
        print ' starting twist: %s  ending twist: %s  number of samples: %s \n'%(w0_begin,w0_end,w0_iter)
        print ' starting R0:   %s   ending R0:  %s   number of smaples:  %s \n'%(R0_begin,R0_end,R0_iter)
        print ' number of chains:  %s \n'%(num_chain)
        print ' number of chains in assymetric unit: %s \n' %(num_to_output)
        print ' ###############  individual helix parameters  ######## \n'
        print ' ORIENTATION CHAIN CHAIN_ORDER   PHASE (start, end, number)   Z OFFSET (start, end, number) \n '
        for i in range(num_to_output):
            print '%s %s %s    %s %s %s              %s %s %s \n'%(orientation[i],chain_name[i],chain_order[i],phase[i][0],phase[i][1],phase[i][2],z_list[i][0],z_list[i][1],z_list[i][2])

        number_of_combinations=w0_iter*R0_iter
        for i in range(num_to_output):
            number_of_combinations=number_of_combinations*phase[i][2]*z_list[i][2]

        print  '  NUMBER OF PDBS TO GENERATE:  %s \n'%number_of_combinations


        combinations=[]
        w0_inc = (w0_end - w0_begin)/max(w0_iter-1,1)
        R0_inc = (R0_end - R0_begin)/max(R0_iter-1,1)
        w0=[]
        R0=[]
        ph=[]
        z = []
        for i in range(w0_iter):
            w0.append( w0_begin+w0_inc*i)
            # insert computation of w1
            # PG

        for i in range(R0_iter):
            R0.append( R0_begin+R0_inc*i)

        for j in range(num_to_output):
            p_inc=(phase[j][1] - phase[j][0])/max(phase[j][2]-1,1)
            ph.append([])
            for i in range(phase[j][2]):
                ph[j].append(phase[j][0]+p_inc*i)

        for j in range(num_to_output):
            z_inc=(z_list[j][1] - z_list[j][0])/max(z_list[j][2]-1,1)
            z.append([])
            for i in range(z_list[j][2]):
                z[j].append(z_list[j][0]+z_inc*i)

        return(Nres,num_chain,num_to_output,orientation,R0,w0,ph,z,output_file_prefix,chain_name,chain_order)

####################################################                      

    def main(self):

        input_file=argv[1]
        tag = argv[2]

        input =map(string.split,open(input_file,'r').readlines())
        Nres,num_chain,num_to_output,orientation,R0,w0,ph,z,output_file_prefix,chain_name,chain_order = self.input_params(input)
        ##print w0,ph,z,chain_order

        items=[]
        items.append(len(w0))
        items.append(len(R0))

        for p in ph:
            items.append(len(p))
        for zz in z:
            items.append(len(zz))

        combos=self.enumerate_combinations(items)

        for combo in combos:
            id=map(int,string.split(combo))
            helix_phase=[]
            delta_z =[]

            for i in range(num_to_output):
                helix_phase.append(ph[i][id[2+i]])
                delta_z.append(z[i][id[2+num_to_output+i]])
    
            out_file_name='%s_%.2f_%.2f'%(tag,w0[id[0]],R0[id[1]])
            for i in range(num_to_output):
                out_file_name=out_file_name+'_'+'%.2f'%helix_phase[i]
            for i in range(num_to_output):
                out_file_name=out_file_name+'_'+'%.2f'%delta_z[i]
            out_file_name=out_file_name+'.pdb'
            self.Generate_Pdbs(Nres,num_chain,num_to_output,w0[id[0]],R0[id[1]],orientation,helix_phase,delta_z,out_file_name,chain_name,chain_order)


if __name__ == "__main__":
    run = ParametricHelicalBundle()
    run.main()

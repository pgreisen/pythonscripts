#!/usr/bin/env python
from os import system,popen
import string
from sys import argv
import math
from math import sin,cos
#import xyzMath
from xyzMath import *


def enumerate_combinations(list):
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

def assign_heptad(ph):
    std=[197,300,41,146,249,351,95]
    label=['a','b','c','d','e','f','g']
    while ph < 0. :
        ph=ph+360.
    while ph > 360. :
        ph=ph-360.
    d_c=200.**2
    close=8
    for n in range(7):
       d=(ph-std[n])**2
       if d < d_c:
           d_c =d
           close=n
   # print ph, close
    return(label[close])

def Generate_Pdbs(Nres,num_chain,num_to_output,w,R,orientation,helix_phase,delta_z,output_file_name,chain_name, chain_order):

 deg_to_rad= math.pi/180.
# chain parameters
 ph = 360/num_chain
 phase=[]
 for i in range(num_chain):
    phase.append(i*ph)
#    phase.append(0)
 chain_set=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

 chain_num=chain_set[0:num_chain]

 R1=2.26
# w1 = 720./7.  # 102.85  for a 2-layer heptad repeat
# w1 = 1080/11
 w1=100.
#rise per residue d  fixed. this constrains pitch (alpha)
 d=1.51
 z1=0.
 z2=0.
 line='ATOM      7  CA  GLY A   2       5.520   2.352   1.361  1.00 20.00     '
 last=line[54:-1]


 atom_num=1
 res_num=1
 Res_id=[]
 CA_list=[]
 heptad=[]
 alpha=math.asin(R*w*deg_to_rad/d)
 for iter in range(num_to_output):
#  iter=int(chain_order[counter])
  CA_list.append([])
  Res_id.append([])
  heptad.append([])
#  chain=chain_num[iter]
  chain = chain_name[iter]
  orient = orientation[iter]

  if orient == 1:
     res_num=0
  else:
     res_num=Nres+1

  supercoil_phase=phase[iter]+delta_z[iter]*math.tan(alpha)/(R*deg_to_rad)
  for t in range(Nres+2):      ## need two extra residues to guide placement of 1st and last residue
    a0=(w*t+supercoil_phase)*deg_to_rad
    #a1=(w1*t+helix_phase[iter]+phase[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius
    a1=(w1*t+helix_phase[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius

    x=R*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
    y=R*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
    if w == 0:
        z=d*t+delta_z[iter]
    else:
        z= R*w*t*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+delta_z[iter]


#    pdb_lines.append( (res_num,'ATOM   %4d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,x,y,z,last)) )

#    print iter
    ph=(w1*t+helix_phase[iter])
    heptad[iter].append( ( res_num,assign_heptad(ph)))
    CA_list[iter].append( (res_num,Vec(x,y,z)) )
    Res_id[iter].append(res_num)
    atom_num=atom_num+1

    if orient  == 1:
        res_num=res_num+1
    else:
        res_num=res_num-1
# convert CA trace to full backbone model by superimposing on ideal template
# by matching 3 consecutive CA atoms
# set up ideal template
 stub_file=map(string.split,open('ideal_ala.pdb','r').readlines())
 atom=[]
 for line in stub_file:
    atom.append( (Vec(float(line[6]),float(line[7]),float(line[8]))))

 ideal_stub=stub(atom[6],atom[1],atom[11])
##  print atom[6]
##  print atom[8]
##  print Vec.distance(atom[6],atom[8]),'dist'
#now make full backbone pdb

 full_pdb=open(output_file_name,'w')
 atom_num=1

 res_n=0
 for counter in range(num_to_output):
  iter=int(chain_order[counter])
  chain=chain_name[iter]
  CA_chain_u=CA_list[iter]
  heptad_chain_u=heptad[iter]
  heptad_chain=sorted(heptad_chain_u, key = lambda res: res[0])
  CA_chain = sorted(CA_chain_u, key = lambda res: res[0])
  for res in range(1,Nres+1):
    res_n=res_n+1
    #print res_n,heptad_chain[res][1]
    res_num='%3d%s'%(res_n,heptad_chain[res][1])
    #print res_num
#    print res, res_num
#    print res_num
#    print CA_chain[res][1]
    actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])
    transform=actual_stub * ~ideal_stub

#    start_d=Vec.distance(atom[5],atom[6])
#    end_d = Vec.distance(coords,transform*atom[6])
#    print start_d,end_d,'dist'
#    print CA_list[res],'ori',coords

# N
    coords=transform*atom[5]
    full_pdb.write('ATOM %6d  N   ALA %s%s    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

# CA   (use actual CA from trace rather than superimposed one)
    coords=CA_chain[res][1]
    tcoords=transform*atom[6]
#    print coords,tcoords,'CA'

    full_pdb.write('ATOM %6d  CA  ALA %s%s    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

#  C
    coords=transform*atom[7]
    full_pdb.write('ATOM %6d  C   ALA %s%s    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

# O
    coords=transform*atom[8]
    full_pdb.write('ATOM %6d  O   ALA %s%s    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

    start_d=Vec.distance(atom[8],atom[6])
    end_d = Vec.distance(transform*atom[8],transform*atom[6])
# CB
    coords=transform*atom[9]
    full_pdb.write('ATOM %6d  CB  ALA %s%s    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

#    print start_d,end_d,'dist C  CA '
 return()


def input_params(input):

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
    # sample p1 evenly, w and R around the input values

    ## z_list=[1.5,2.0,2.5,3.0]
    ## p1_s=125.
    ## p2_s=235.


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

input_file=argv[1]
tag = argv[2]
input =map(string.split,open(input_file,'r').readlines())
Nres,num_chain,num_to_output,orientation,R0,w0,ph,z,output_file_prefix,chain_name,chain_order = input_params(input)
#print w0,ph,z,chain_order
items=[]
items.append(len(w0))
items.append(len(R0))
for p in ph:
    items.append(len(p))
for zz in z:
    items.append(len(zz))

combos=enumerate_combinations(items)

for combo in combos:
#    print combo
    id=map(int,string.split(combo))
#    print id
    helix_phase=[]
    delta_z =[]

    for i in range(num_to_output):
        helix_phase.append(ph[i][id[2+i]])
        delta_z.append(z[i][id[2+num_to_output+i]])
    ## to get close to C2 symmetry for axis in x-y plane, switch sign of  phase and delta_z of sym mates
        ##
##     if num_chain==2:
## 	helix_phase.append(orientation[1]*helix_phase[0])
##         delta_z.append(orientation[1]*delta_z[0])

##     if num_chain==4:   ## warning-this is 4 helix bundle centric
## 	helix_phase.append(-helix_phase[1])
##         helix_phase.append(-helix_phase[0])
##         delta_z.append(-delta_z[1])
##         delta_z.append(-delta_z[0])




    out_file_name='%s_%.2f_%.2f'%(tag,w0[id[0]],R0[id[1]])
    #out_file_name='%s_%.2f'%(tag,w0[id[0]])
    for i in range(num_to_output):
        out_file_name=out_file_name+'_'+'%.2f'%helix_phase[i]
    #for i in range(num_to_output):
    #    out_file_name=out_file_name+'_'+'%.2f'%delta_z[i]
    out_file_name=out_file_name+'.pdb'
    Generate_Pdbs(Nres,num_chain,num_to_output,w0[id[0]],R0[id[1]],orientation,helix_phase,delta_z,out_file_name,chain_name,chain_order)

## for i in range(3):
##     p1=p1_s+(i-1)*5
##     for ii in range(3):
##      p2=p2_s+(ii-1)*5
##      for j in range(3):
##         R=R0+(j-1)*0.1
##         for k in range(5):
##             w=w0+(k-2)*0.15
##             for l in range(4):
##                 delta_z=z_list[l]
##                 output_file='fine2_%s_%s_%s_%s_%s_%s'%(w,R,p1,p2,delta_z,output_file_name)
## #    print output_file_name,Nres,w,R,p1
##                 Generate_Pdbs(Nres,w,R,num_chain,p1,p2,delta_z,output_file,num_to_output)
## #print atom_list

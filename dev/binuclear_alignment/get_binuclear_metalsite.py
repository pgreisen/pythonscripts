from numpy import *
from math import sqrt
import os,shutil

class get_binuclear_metalsite:

    # Requires PDB files
    # Returns object of file
    def get_obj_pdbfile(self,pdbfile):
        fl = open(pdbfile,'r')
        pdb = fl.readlines()
        fl.close()
        return pdb

    # Requires name of metal as string
    # Returns a dictionary with the coordinates
    def get_metal_ions(self,pdbfile,mt=['ZN','MG','MN','FE','CU','NI','CO','CA','CD']):
        ht_lig = ['O','S','N']
        coor_zn = {}
        hetatm = {}
        for line in pdbfile:
            if line[0:6] == 'HETATM':
                tmp = line.split()
                het_id = line[13:14]

                # We need to include an extra verification parameter due to
                # CA as calcium but also CA for other ligands so we try the CA for atom
                # and CA for residue

                metal_resname = str(line[17:20]).strip()
                if tmp[2] in mt and metal_resname in mt:
                    tmp_id = line[22:26] 
                    coor_zn[tmp_id] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                elif  het_id in ht_lig:
                    nid = line[7:11]+'   '+line[13:16]   #het_id
                    hetatm[nid] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
        return coor_zn,hetatm
    

    # Getting the distance between zinc ion and hetero atoms
    # returning the atoms less then DISTANCE from the ZN as well
    # as a boolean value in order to search for protein ligands
    # Requires dictionary with metal coordinates
    # Returns dictionary with active site as well as boolean
    # 04-11-2009 remove boolean as geometry is correct
    def get_hetero_around_metal(self,metal_dic,het_dic,DISTANCE):
        active_site = {}
        for i in metal_dic:
            metal_vec = metal_dic[i]
            for j in het_dic:
                value =  het_dic[j]
                metaldis = linalg.norm(metal_vec - value)
                if metaldis <= DISTANCE:
                    active_site[j]=value
        return active_site


    # getting single metal site
    # metal_dic is tuple
    def get_hetero_around_single_metal(self,metal_dic,het_dic):
        active_site = {}
        # metal_dic with name and vector
        metal_vec = metal_dic[1]
        for j in het_dic:
            value =  het_dic[j]
            metaldis = linalg.norm(metal_vec - value)
            if metaldis < 2.7:
                active_site[j]=value
        return active_site



    

    # Requires fileobject, metal coordinates, and active site coordinates
    # Returns dictionary with ligands <= DISTANCE  from Zn,list of protein ligands
    # and their distance to Zn
    def get_protein_ligand_metal(self,pdb_obj,zn_coor,DISTANCE): #,active_site):
        # Distance between ligands and metal ions
        mt_lig = []
        active_site = {}
        coor_residues = ['CYS','HIS','ASP','GLU','SER','GLN']
        lig_atm = ['ND1', 'NE2','SG ', 'OE1', 'OE2','OD1','OD2','OG ']
        for line in pdb_obj:
            res = line[17:20]
            atm = line[0:4]
            if res in coor_residues and atm =='ATOM':
                atom = line[13:16]
                if atom in lig_atm:
                    
                    vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                    for i in zn_coor:
                        ds_lig = linalg.norm(vec-zn_coor[i])
                        if linalg.norm(vec-zn_coor[i]) <= DISTANCE:
                            iid = res+'  '+line[23:26]
                            mt_lig.append(iid+'  '+atom+'  '+str(ds_lig))
                            active_site[iid] = vec
        return active_site


    # Requires fileobject, metal coordinates, and active site coordinates
    # Returns dictionary with ligands < 2.5 from Zn,list of protein ligands
    # and their distance to Zn
    def get_protein_ligand_single_metal(self,pdb_obj,zn_coor,active_site,DISTANCE):
        # Distance between ligands and metal ions
        mt_lig = []
        coor_residues = ['CYS','HIS','ASP','GLU','SER','GLN','ASN']
        lig_atm = ['ND1', 'NE2','SG ', 'OE1', 'OE2','OD1','OD2','OG ']
        for line in pdb_obj:
            res = line[17:20]
            atm = line[0:4]
            if res in coor_residues and atm =='ATOM':
                atom = line[13:16]
                if atom in lig_atm:
                    vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                    ds_lig = linalg.norm(vec-zn_coor)
                    if ds_lig < DISTANCE:
                        iid = res+'  '+line[23:26]
                        mt_lig.append(iid+'  '+atom+'  '+str(ds_lig))
                        active_site[iid] = vec
        return active_site,mt_lig

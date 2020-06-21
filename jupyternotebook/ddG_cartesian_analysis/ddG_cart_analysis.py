#!/usr/bin/env python
# coding: utf-8

# # Analysis of ddG computation from Rosetta
xtal_off = 0


import os
import pandas as pd
import glob
import matplotlib.pyplot as plt
from scipy import stats as scistats
import numpy as np


# get fasta file from input pdb-structure used for ddG calculations to ensure 
# this will make assert on the merge with PSSM and ddG more secure


def get_single_letter_AA(residue):
    aa = {
        "ALA" : "A",
        "ILE" : "I",
        "LEU" : "L",
        "VAL" : "V",
        "MET" : "M",
        "PHE" : "F",
        "TYR" : "Y",
        "ARG" : "R",
        "LYS" : "K",
        "TRP" : "W",
        "ASN" : "N",
        "GLN" : "Q",
        "ASP" : "D",
        "GLU" : "E",
        "SER" : "S",
        "THR" : "T",
        "GLY" : "G",
        "PRO" : "P",
        "CYS" : "C",
        "HIS" : "H"
        }
    if(isinstance(residue,str)):
        return aa[residue]
    return "NaN"


def get_stats(group):
    return {'std': np.std(group),'mean': group.mean()}

def get_ddG_monomer_results(path_and_file_name):
    header_text=[
    'COMPLEX:', 'Round', 'Baseline', 'Total Energy', 'fa_atr_label:', 'fa_atr:',
       'fa_rep_label:', 'fa_rep:', 'fa_sol_label:', 'fa_sol:', 'fa_intra_rep_label:', 'fa_intra_rep:',
       'fa_intra_sol_xover4_label:', 'fa_intra_sol_xover4', 'lk_ball_wtd_label:', 'lk_ball_wtd:', 'fa_elec_label:',
       'fa_elec:', 'hbond_sr_bb_label:', 'hbond_sr_bb:', 'hbond_lr_bb_label:', 'hbond_lr_bb:',
       'hbond_bb_sc_label:', 'hbond_bb_sc:', 'hbond_sc_label:', 'hbond_sc:', 'dslf_fa13_label:',
       'dslf_fa13:', 'omega_label:', 'omega:', 'fa_dun_label:', 'fa_dun:', 'p_aa_pp_label:', 'p_aa_pp:',
       'yhh_planarity_label:', 'yhh_planarity:', 'ref_label:', 'ref:', 'rama_prepro_label:', 'rama_prepro:',
       'cart_bonded_label:', 'cart_bonded:'] 
    df = pd.read_csv(path_and_file_name,skiprows=0,delim_whitespace=True,names=header_text)
    return df



# gather data from computation ordinary ddG computation
dfs = []
directory_w_subdirs = './'
filename="mutfile.ddg"
for dr in glob.glob(directory_w_subdirs+"*"):
    if(os.path.isdir(dr) and dr != 'min_cart'):
        try:
            dfs.append(get_ddG_monomer_results(dr+"/"+filename))
        except:
            print("No ddg file present in directory: "+dr)
df_tot = pd.concat(dfs,sort=False)


# In[10]:


df_tot = df_tot.groupby(["Baseline"])["Total Energy"].apply(get_stats).unstack().reset_index()


# In[11]:


df_tot.rename(columns={"mean" : "total (mean n=10)", 
                       "std" : "std( total )"},inplace=True)


# In[12]:


df_tot["ddG"] = df_tot["total (mean n=10)"] - df_tot["total (mean n=10)"].loc[(df_tot["Baseline"] == "WT:") ].values[0]


df_tot.columns


# In[14]:


#df_tot["Baseline"].unique()
cutoff = df_tot[df_tot["Baseline"] == "WT:"]['ddG'].values[0]+2*df_tot[df_tot["Baseline"] == "WT:"]['std( total )'].values[0]
print(cutoff)


# In[15]:


# position and AA
df_tot["Pos"] = df_tot["Baseline"].str.split("_").str[1].str[0:-4]
df_tot["AA"] = df_tot["Baseline"].str.split("_").str[1].str[-4:-1].apply(get_single_letter_AA)


# In[16]:


df_tot.info()


# In[17]:


plt.figure();
df_tot["ddG"].plot.hist(bins=20)
plt.title("ddG")
plt.xlabel("ddG (kcal/mol)")


# In[18]:


plt.figure()
pd.value_counts(df_tot[df_tot["ddG"] < cutoff]["AA"]).plot.bar()
plt.title("Distributions of substitutions < "+str(round(cutoff,2))+" kcal/mol")


# In[19]:


df_tot.rename(columns={"ddG" : "ddG_cart"},inplace=True)


# In[20]:


df_tot.head()


# In[21]:


df_tot = pd.merge(right=df_wt_seq,right_on="Pos",left=df_tot,left_on="Pos")
df_tot.head()


# In[22]:


df_tot.info()


# In[23]:


print(len(df_tot[df_tot["ddG_cart"] < cutoff]["Pos"].unique()))


# In[24]:


df_tot["AAPos"] = df_tot["Pos"]+df_tot["AA"]
columns_to_file = ["WT_AA","AAPos","Pos","AA","ddG_cart"]
df_tot[columns_to_file].to_csv("20181211_ddG_cart_rosetta_analysis.csv")


# In[32]:


# df_tot.sort_values(by="ddG_cart",ascending=True)

# what is this?
ddg_calc_to_mut = {}
for files in glob.glob("*"):
    for j in df_tot["AAPos"]:
        try:
            if(str(files).endswith(j)):
                ddg_calc_to_mut[j] = files
        except:
            continue
# In[35]:


df_tot.head(1)


# In[47]:


# write to pymol pml-file
with open("ddG_positions.pml",'w') as f:
    f.write("show cartoon\n")
    f.write("hide lines\n")
    dftmp = df_tot[df_tot["ddG_cart"] < cutoff]
    for i,j,k in zip(dftmp["WT_AA"],dftmp["Pos"],dftmp["AA"]):
        #try:
        tmppair = i+str(j)+k
        pair="create "+tmppair
        f.write(pair+",resi "+str(int(j)+xtal_off)+"\n" )
        f.write("show sticks, "+str( tmppair )+"\n")
        f.write("color cyan, "+str( tmppair )+" and name C*\n")
        #except:
        #    print("NaN for wt")
    f.write("show lines\nhide everything, elem h")


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt
import math
import datetime
date = str(datetime.date.today())
from datetime import datetime as dt

import numpy as np
import pandas as pd
from ipywidgets import widgets


# In[2]:


from _plotly_future_ import v4_subplots
from plotly.subplots import make_subplots
import plotly
from plotly.graph_objs import Scatter, Layout
import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff
plotly.offline.init_notebook_mode(connected=False)


# In[3]:


# Inputfile name                                                                                      
inp_file = "2020-01-02_20191219_APassay_protease_challenge_AP-HT_variantsrepeat_project_KineticMeasurements.csv"

# names for output file                                                                               
exp="elanco_20191219_APassay_protease_challenge_AP-HT_variantsrepeat"
read = 2

run="run_"+str(read)+"_"
pth = "/Users/pgreisen/Projects/Enevolv/elanco/20200102_AP_R3R4_HT/"

df_reformat = pd.read_csv(pth+inp_file)
df_reformat = df_reformat[df_reformat['read'] == read]

cols_to_keep = ['temperature', 'well', 'fluorescence', 'time_s']
df_reformat = df_reformat[cols_to_keep]


# In[4]:


df_reformat['fluorescence']  = df_reformat['fluorescence'].replace('OVRFLW',0)
df_reformat['fluorescence'] = df_reformat['fluorescence'].astype(float)


# In[5]:


df = df_reformat.pivot_table(index=['time_s', 'temperature'],columns=['well'])
df.columns = df.columns.droplevel(0)
df = df.rename_axis(None, axis=1)
df = df.reset_index()
df = df.drop(labels=['temperature'],axis=1)
order_of_output_columns = sorted(list(df.columns), key=lambda x: (x[0], x[1:-1]))
# this has been changed to seconds and is converted into minutes here                                 
df["time_min"] = df['time_s'] / 60
df.fillna(value=0.0,inplace=True)


# In[6]:


dataframes = df

yvalue = widgets.Dropdown(
    description='WellID:   ',
    value=list(df.columns)[0],
    options=list(df.columns)
)


# In[7]:


xvalue = widgets.Dropdown(
    description='Time:   ',
    value='time_s',
    options=['time_s','time_min'] )


# In[8]:


container = widgets.HBox(children=[xvalue,yvalue])


# In[9]:


trace_scatter = go.Scatter(
    x=dataframes[xvalue.value],
    y=dataframes[yvalue.value],
    mode='lines+markers')


# In[10]:


g_scatter = go.FigureWidget(data=[trace_scatter],
                    layout=Layout(title=yvalue.value,
                           xaxis=dict(title=xvalue.value,
                                                  titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f')),
                           yaxis=dict(
                               title=yvalue.value,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f'))))


# In[11]:


def validate():
    if yvalue.value in list(dataframes.keys()):
        return True
    else:
        return False


# In[12]:


def response(change):
    if validate():
        x1 = dataframes[xvalue.value]
        y1 = dataframes[yvalue.value]
        with g_scatter.batch_update():
            g_scatter.data[0].x = x1
            g_scatter.data[0].y = y1
            g_scatter.layout.xaxis.title = xvalue.value
            g_scatter.layout.yaxis.title = yvalue.value


# In[13]:


xvalue.observe(response, names="value")
yvalue.observe(response, names="value")
widgets.VBox([container,g_scatter])


# In[ ]:





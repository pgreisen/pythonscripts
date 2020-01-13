#!/usr/bin/env python
# coding: utf-8

# ## Example of jupyterlab-dash extension

# In[5]:


# Imports
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go


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
#from _plotly_future_ import v4_subplots
from plotly.subplots import make_subplots
import plotly
from plotly.graph_objs import Scatter, Layout
import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff


# In[7]:


# Build AppViewer 
from jupyterlab_dash import AppViewer
viewer = AppViewer()


# In[8]:


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
dataframes = df


# In[4]:
import sys                                                                                           
sys.path.append("/Users/pgreisen/pythonscripts/develop/Dashboards/LinearFit")                        
import LinearFit
lf = LinearFit.LinearFit()
lf.main()
tmp = lf.get_fitted_data()

first_key = list(tmp.keys())[]





# Build App
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.Div([
        html.Div([
            dcc.Dropdown(
                id='xaxis',
                options=[{'label': i, 'value': i} for i in ['time_s','time_min']],
                value='time_s'
            )
        ],
        style={'width': '48%', 'display': 'inline-block'}),
        html.Div([
            dcc.Dropdown(
                id='yaxis',
                options=[{'label': i, 'value': i} for i in list(df.columns)],
                value='A1'
            ),],style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),

    dcc.Graph(id='indicator-graphic'),
    # Adding more to the dashboard
    html.Div([
            dcc.Dropdown(
                id='xaxis2',
                options=[{'label': i, 'value': i} for i in list(df.columns)],
                value='A1'
            ),],style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

    html.Div([
            dcc.Dropdown(
                id='yaxis2',
                options=[{'label': i, 'value': i} for i in list(df.columns)],
                value='A1'
            ),],style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),


    
    ])

# Callbacks
@app.callback(
    dash.dependencies.Output('indicator-graphic', 'figure'),
    [dash.dependencies.Input('xaxis', 'value'),
     dash.dependencies.Input('yaxis', 'value')])

def update_graph(xaxis_column_name, yaxis_column_name):
    dff = df
    return {
        'data': [go.Scatter(
            x=dff[xaxis_column_name],
            y=dff[yaxis_column_name],
            mode='lines+markers',
        )],
        'layout': go.Layout(
            xaxis={
                'title': xaxis_column_name},
            yaxis={
                'title': yaxis_column_name,},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
            hovermode='closest',
        )
    }

# to view inside the lab
#viewer.show(app)


# In[21]:


if __name__ == '__main__':
    app.run_server(debug=True)


# In[ ]:





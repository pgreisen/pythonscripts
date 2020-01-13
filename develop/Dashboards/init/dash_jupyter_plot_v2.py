#!/usr/bin/env python
# coding: utf-8

# ## Example of jupyterlab-dash extension

# In[5]:


# Imports
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
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

# Build AppViewer 
from jupyterlab_dash import AppViewer
viewer = AppViewer()

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
dfs2, fitted_data = lf.get_fitted_data()
keys_2 = list(dfs2.keys())

dfout = pd.read_excel( lf.outputname, sheet_name=None )

values2 = ['alpha','b','Residual','R2','Rank','CV']

# Build App
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# Dash app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


app.layout = html.Div([
    html.Div([
        html.H5("Raw data"),
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
        html.H5("Scatter plot of data"),

        html.Div([
            dcc.Dropdown(
                id='exp',
                options=[{'label': i, 'value': i} for i in keys_2],
                value=keys_2[0]
            ), ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

        html.Div([
            dcc.Dropdown(
                id='exp_fit',
                options=[{'label': i, 'value': i} for i in values2],
                value=values2[0]
            ), ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

    html.Div([
        dcc.Dropdown(
                id='exp_fit2',
                options=[{'label': i, 'value': i} for i in values2],
                value=values2[0]
            ),],style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),


        ]),
html.Br([]),html.Br([]),html.Br([]),

    html.Div([
        dcc.Graph(id='graphB'), ]),

    html.Br([]),html.Br([]),html.Br([]),

    html.Div([
html.H5("Linear fit of data"),

        html.Div([
            dcc.Dropdown(
                id='name_dropdown',
                options=[{'label': i, 'value': i} for i in list(fitted_data.keys())],
                value=list(fitted_data.keys())[0]
            ), ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

        html.Div([
            dcc.Dropdown(
                id='opt_dropdown',), ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),
        html.Hr(),
        html.Div(id='display-selected-values'),]), # End section H5
    html.Div([
        dcc.Graph(id='graphC'), ]),

html.Br([]), html.Br([]), html.Br([]),
    html.Br([]), html.Br([]), html.Br([]),

    html.Div([
    html.H5("Slope of fits"),
        dcc.Dropdown(
            id='table_display_exp',
            options=[{'label': i, 'value': i} for i in list(dfout.keys())],
            value='Summary'
        ), ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

    html.Br([]), html.Br([]), html.Br([]),

    html.Div([
        dash_table.DataTable(id='table',
                             columns=[{"name": i, "id": i} for i in dfout['Summary'].columns],
                             data=dfout['Summary'].to_dict('records'),
                             ), ]
    ),



    html.Div(
    html.H5("Top slope for well"),
        id='output-data-upload'),
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
                'title': yaxis_column_name, },
            margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
            hovermode='closest',
        )
    }


# Callbacks
@app.callback(
    dash.dependencies.Output('graphB', 'figure'),
    [dash.dependencies.Input('exp', 'value'),
     dash.dependencies.Input('exp_fit', 'value'),
     dash.dependencies.Input('exp_fit2', 'value')])


def update_graphB(exp, exp_fit, exp_fit2):
    dff = dfs2
    return {
        'data': [go.Scatter(
            x=dff[exp][exp_fit],
            y=dff[exp][exp_fit2],
            mode='markers',
            marker=dict(
            size=16,
            color=dff[exp]['R2'],  # set color equal to a variable
            colorscale='Viridis',  # one of plotly colorscales
            showscale=True))],
        'layout': go.Layout(
            xaxis={
                'title': exp_fit},
            yaxis={
                'title': exp_fit2, },
            margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
            hovermode='closest',
        )
    }



@app.callback(
    dash.dependencies.Output('opt_dropdown', 'options'),
    [dash.dependencies.Input('name_dropdown', 'value')]
)
def update_date_dropdown(name):
    return [{'label': i, 'value': i} for i in fitted_data[name]]

@app.callback(
    dash.dependencies.Output('display-selected-values', 'children'),
    [dash.dependencies.Input('opt_dropdown', 'value')])

def set_display_children(selected_value):
    return 'you have selected {} option'.format(selected_value)


@app.callback(
    dash.dependencies.Output('graphC', 'figure'),
    [dash.dependencies.Input('name_dropdown', 'value'),
     dash.dependencies.Input('opt_dropdown', 'value')])


def update_graphC(name_dropdown, opt_dropdown):
    dff = fitted_data
    x1 = dff[name_dropdown][opt_dropdown][-4]
    y1 = dff[name_dropdown][opt_dropdown][-3]
    y2 = dff[name_dropdown][opt_dropdown][-2]
    trace1 = go.Scatter(
        x=x1,
        y=y1,
        mode='lines', name='Raw data')

    trace2 = go.Scatter(
        x=x1,
        y=y2,
        mode='lines',name='Fitted line')

    return {'data': [trace1, trace2],'layout': go.Layout(xaxis={'title': 'test'},
            yaxis={
                'title': 'test', },
            margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
            hovermode='closest',
        )
    }


@app.callback(dash.dependencies.Output('output-data-upload', 'children'),
			 [dash.dependencies.Input('table_display_exp', 'value'),])


def update_output(table_display_exp):
	return html.Div([
        dash_table.DataTable(
		id='table2',
        columns=[{"name": i, "id": i} for i in dfout[table_display_exp].columns],
        data=dfout[table_display_exp].to_dict('records'))])

'''
@app.callback(
    dash.dependencies.Output('table','table'),
    [dash.dependencies.Input('table_display_exp', 'value')])


def update_table(table_display_exp):
    return dfout[table_display_exp]
'''
'''
@app.callback(dash.dependencies.Output('table2', 'children'),
			 [dash.dependencies.Input('table_display_exp', 'value')])

def update_output( table_display_exp):
    return dfout[table_display_exp]
'''

# to view inside the lab
#viewer.show(app)


# In[21]:


if __name__ == '__main__':
    app.run_server(debug=True)


# In[ ]:





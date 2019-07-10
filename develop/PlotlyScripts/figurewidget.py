import pandas as pd
import plotly.graph_objs as go
import ipywidgets as widgets

items = df.columns.unique().tolist()
x = widgets.Dropdown(options=items, description='X-value')
y = widgets.Dropdown(options=items, description='Y-value')


def response(change):
    xtmp = x.value
    ytmp = y.value
    
    x_ = df[xtmp]
    y_ = df[ytmp]
    
    fig.data[0].x = x_
    fig.data[0].y = y_
    
# Assign an emptry figure widget with two traces
trace1 = go.Scatter(x=df['sc'],y=df['LigIntFilt'], opacity=0.75, name='Rosetta Values',\
                   mode='markers',text=df['description'], marker=dict(size=10))

fig = go.FigureWidget(data=[trace1])  
x.observe(response, names='value')
display(widgets.VBox([x,y, fig]))

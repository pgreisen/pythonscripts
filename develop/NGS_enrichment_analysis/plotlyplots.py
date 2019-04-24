import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff

def barplot(df,cols,title):
    data = []
    for i in cols:
        data.append( go.Bar(x=[i],y=df[i],name=str(i) ) )
    layout = go.Layout(title=title)
    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig, filename="test.html",auto_open=False)

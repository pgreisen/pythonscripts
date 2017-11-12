import plotly
from plotly.graph_objs import Scatter, Layout

plotly.offline.init_notebook_mode(connected=False)

plotly.offline.iplot({
    "data": [Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])],
    "layout": Layout(title="hello world")
})

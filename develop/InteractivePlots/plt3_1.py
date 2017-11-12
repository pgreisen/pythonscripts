import plotly
from plotly.graph_objs import Scatter, Layout


import plotly.graph_objs as go
import plotly.plotly as py
import plotly
import numpy as np

import plotly
#plotly.offline.init_notebook_mode() # run at the start of every ipython notebook
#plotly.offline.iplot()


trace1 = go.Scatter(
    y = np.random.randn(500),
    mode='markers',
    marker=dict(
        size='16',
        color = np.random.randn(500), #set color equal to a variable
        colorscale='Viridis',
        showscale=True
    )
)

data = [trace1]

plotly.offline.plot({
    "data": data,
    "layout": Layout(title="hello world")
})

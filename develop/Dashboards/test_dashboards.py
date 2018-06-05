import plotly.dashboard_objs as dashboard
import IPython.display
from IPython.display import Image
 
my_dboard = dashboard.Dashboard()
my_dboard.get_preview()


box_1 = {
    'type': 'box',
    'boxType': 'plot',
    'fileId': 'PlotBot:1296',
    'title': 'scatter-for-dashboard'
}
 
box_2 = {
    'type': 'box',
    'boxType': 'plot',
    'fileId': 'PlotBot:1298',
    'title': 'pie-for-dashboard'
}
 
box_3 = {
    'type': 'box',
    'boxType': 'plot',
    'fileId': 'PlotBot:1342',
    'title': 'box-for-dashboard',
    'shareKey':'uVoeGB1i8Xf4fk0b57sT2l'
}
 
my_dboard.insert(box_1)
my_dboard.insert(box_2, 'above', 1)
my_dboard.insert(box_3, 'left', 2)

import plotly.plotly as py
py.dashboard_ops.upload(my_dboard, 'My First Dashboard with Python')

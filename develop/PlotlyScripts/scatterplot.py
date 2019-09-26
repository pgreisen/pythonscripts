def get_scatter_plot(dataframe,xlabel,ylabel,description='description',plottitle="Default"):
    data = []
    data.append(go.Scatter(x=dataframe[xlabel],y=dataframe[ylabel],mode='markers',\
                           text=dataframe[description], marker=dict(size=10)))
    layout = go.Layout(title=plottitle,xaxis=dict(title=xlabel,
                                                  titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f')),
                           
                                                     yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))
    fig = go.Figure(data=data, layout=layout)
    return fig

var1= "NAME"
fig = get_scatter_plot(df, x, y,'description',var1 )
offline.iplot(fig, filename=var1+".png",image='png')
offline.plot(fig, image='png', filename=var1+".html", image_filename=var1,auto_open=False)

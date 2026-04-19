import numpy as np
import plotly.graph_objects as go

def plot_2d(y:list,x:list = None, xlabel:str = None, ylabel:str = None,title:str = None,line_name:list[str] = None)-> None:
    """ 
    Do a simple 2D line chart (list of lines or other)
    """
    fig = go.Figure()
    

    if type(y[0]) == list or type(y[0]) == np.ndarray:
        for i in range(len(y)):
            if x is not None:
                fig.add_trace(go.Scatter(
                    x=x[i],
                    y=y[i],
                    mode='lines',
                    name=line_name[i]
                ))
            else:
                fig.add_trace(go.Scatter(
                    y=y[i],
                    mode='lines',
                    name=line_name[i]
                ))
    else:
        if x is not None:
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode='lines',
                name=line_name
            ))
        else:
            fig.add_trace(go.Scatter(
                y=y,
                mode='lines',
                name=line_name
            ))

    if title:
        fig.update_layout(
        title=title)
    if xlabel:
        fig.update_layout(
            xaxis_title = xlabel
        )
    if ylabel:
        fig.update_layout(
            yaxis_title = ylabel
        )

    # 4. Display the figure
    fig.show()
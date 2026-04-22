import numpy as np
import plotly.graph_objects as go
from misc_functions import vector_extract
from attitude import Attitude
from reference_frames import Reference

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

    import plotly.io as pio
    pio.renderers.default = "notebook"
    # 4. Display the figure
    fig.show()

def plot_2d_with_schedule(y:list,x:list = None, xlabel:str = None, ylabel:str = None,title:str = None,line_name:list[str] = None,control_schedule:list[dict]=None,percent=0.001)-> None:
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

    if control_schedule is not None:
        for event in control_schedule:
            fig.add_vline(x=event['t_start'])
            fig.add_vline(x=event['t_end'])
            fig.add_annotation(
        # xref="paper", yref="paper",
        x=event['t_start'] + percent*max(x[0]), y= 0.75,
        textangle=-90,  
        text=f"{event['ref']}",
        showarrow=False,
        font=dict(size=12, color="gray")
    )
        # 4. Display the figure
    import plotly.io as pio
    pio.renderers.default = "notebook"
    fig.show()


def make_full_paper_plots(ref,t0=0,tf=500):
    """ 
    Make all the plots for the paper so they look nice
    """

    # Setup
    lmo = Attitude('lmo')
    X0 = np.hstack((lmo.mrp_bn_0,lmo.omega_bn_0))
    tspan = [t0,tf]

    # Control setup
    P = 1/6
    K = 1/180
    control_schedule = [{'ref':ref,'t_start':t0,'t_end':tf,'K':K,'P':P}]
    lmo.control_schedule = control_schedule

    # Integrating
    attitude_simulated_controlled = lmo.attitude_propagation(X0,tspan,dt=1)

    # Plotting state
    sxc,syc,szc = vector_extract(attitude_simulated_controlled['sigma'])
    tc = attitude_simulated_controlled['time']
    name = [r'$\sigma_1$',r'$\sigma_2$',r'$\sigma_3$']
    plot_2d(y=[sxc,syc,szc],x=[tc,tc,tc],line_name=name,title=f"Simulated MRP - {ref.capitalize()} Pointing",xlabel='Time (s)',ylabel=r'$\sigma_i$')

    # plotting error
    sxc,syc,szc = vector_extract(attitude_simulated_controlled['sigma_error'])
    tc = attitude_simulated_controlled['time']
    plot_2d(y=[sxc,syc,szc],x=[tc,tc,tc],line_name=name,title=f'Simulated MRP Error - {ref.capitalize()} Pointing',xlabel='Time (s)',ylabel=r'$\sigma_i$')

    # Plotting Control
    sxc,syc,szc = vector_extract(attitude_simulated_controlled['u'])
    tc = attitude_simulated_controlled['time']
    plot_2d(y=[sxc,syc,szc],x=[tc,tc,tc],line_name=name,title=f'Simulated Control Vector - {ref.capitalize()} Pointing',xlabel='Time (s)',ylabel=r'$u(X,t)\; (Nm)$')

    # Plotting Omega
    # Plotting Control
    reference = Reference('lmo')
    omega_ref = reference.gmo_pointing_ref(tf,terminal_out=True)
    sxc,syc,szc = vector_extract(attitude_simulated_controlled['omega'])
    tc = attitude_simulated_controlled['time']
    plot_2d(y=[sxc,syc,szc],x=[tc,tc,tc],line_name=name,title=f'Simulated Omega  - {ref.capitalize()} Pointing',xlabel='Time (s)',ylabel=r'$\omega \; (rad/s)$')

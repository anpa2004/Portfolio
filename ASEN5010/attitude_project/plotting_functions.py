import numpy as np
import plotly.graph_objects as go
from misc_functions import vector_extract
from attitude import Attitude
from reference_frames import Reference
from PIL import Image

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


def plot_mars_3d_with_frames(
    trajectories,
    frames=None,
    title="Mars Orbit Visualization",
    mars_radius=3390,
    mars_texture_path=None,
    trajectory_labels=None,
    trajectory_colors=None,
    frame_scale=500,
    show_mars=True,
    show_axes=False,
    show_body_frame_legend=True,
    add_inertial_frame=False,
    inertial_scale=None,
    label_inertial=True
):
    """
    Full 3D Mars visualization with trajectories and coordinate frames.

    Args:
        trajectories : array OR list of arrays (N,3)
        frames : list of dicts
        Each dict:
        {
            "origin": (3,),
            "dcm": (3,3),
            "scale": float (optional),
            "label": str (optional, e.g. 'b'),
            "label_axes": bool (optional)
        }

        add_inertial_frame : bool
            Adds inertial frame at origin

        label_inertial : bool
            Labels inertial frame axes at tips (n̂₁, n̂₂, n̂₃)
    """

    # ------------------------------------------------------------
    # Data Handling
    # ------------------------------------------------------------
    if isinstance(trajectories, np.ndarray):
        trajectories = [trajectories]

    n_traj = len(trajectories)

    if trajectory_labels is None:
        trajectory_labels = [f"Trajectory {i+1}" for i in range(n_traj)]

    if trajectory_colors is None:
        trajectory_colors = [None] * n_traj

    fig = go.Figure()

    # ------------------------------------------------------------
    #  Mars Sphere
    # ------------------------------------------------------------
    if show_mars:
        n = 100
        theta = np.linspace(0, 2*np.pi, n)
        phi = np.linspace(0, np.pi, n)
        theta, phi = np.meshgrid(theta, phi)

        x = mars_radius * np.sin(phi) * np.cos(theta)
        y = mars_radius * np.sin(phi) * np.sin(theta)
        z = mars_radius * np.cos(phi)

        if mars_texture_path:
            img = np.asarray(Image.open(mars_texture_path))
            from skimage.transform import resize
            img = resize(img, (n, n), anti_aliasing=True)

            surfacecolor = np.arange(n*n).reshape(n, n)
            flat_img = img.reshape(-1, 3)

            colorscale = [
                [i/(len(flat_img)-1),
                 f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})']
                for i, (r, g, b) in enumerate(flat_img)
            ]

            fig.add_surface(
                x=x, y=y, z=z,
                surfacecolor=surfacecolor,
                colorscale=colorscale,
                showscale=False,
                hoverinfo='skip'
            )
        else:
            fig.add_surface(
                x=x, y=y, z=z,
                colorscale='Reds',
                showscale=False,
                hoverinfo='skip',
                lighting=dict(ambient=0.7, diffuse=0.8)
            )

    # ------------------------------------------------------------
    #  Trajectories
    # ------------------------------------------------------------
    for i, traj in enumerate(trajectories):
        traj = np.asarray(traj)

        fig.add_trace(go.Scatter3d(
            x=traj[:, 0],
            y=traj[:, 1],
            z=traj[:, 2],
            mode='lines',
            name=trajectory_labels[i],
            line=dict(width=4, color=trajectory_colors[i])
        ))

    # ------------------------------------------------------------
    # Legend (Body Axes)
    # ------------------------------------------------------------
    if show_body_frame_legend:
        axis_info = [
            ("b₁ axis", "red"),
            ("b₂ axis", "green"),
            ("b₃ axis", "blue"),
        ]

        for name, color in axis_info:
            fig.add_trace(go.Scatter3d(
                x=[None], y=[None], z=[None],
                mode='lines',
                line=dict(color=color, width=6),
                name=name,
                showlegend=True
            ))

    # ------------------------------------------------------------
    #  Frames
    # ------------------------------------------------------------
    if frames is None:
        frames = []

    # Add inertial frame if requested
    if add_inertial_frame:
        if inertial_scale is None:
            inertial_scale = mars_radius * 1.5

        frames.append({
            "origin": np.array([0, 0, 0]),
            "dcm": np.eye(3),
            "scale": inertial_scale,
            "label": "n",
            "label_axes": label_inertial
        })

    # Plot frames
    for frame in frames:
        origin = np.asarray(frame["origin"])
        dcm = np.asarray(frame["dcm"])
        scale = frame.get("scale", frame_scale)

        label_prefix = frame.get("label", "b")
        label_axes = frame.get("label_axes", False)

        axes = [
            ("x", "red",   dcm[:, 0], 1),
            ("y", "green", dcm[:, 1], 2),
            ("z", "blue",  dcm[:, 2], 3),
        ]

        for name, color, vec, idx in axes:
            end = origin + scale * vec

            # Draw axis line
            fig.add_trace(go.Scatter3d(
                x=[origin[0], end[0]],
                y=[origin[1], end[1]],
                z=[origin[2], end[2]],
                mode='lines',
                line=dict(color=color, width=6),
                showlegend=False
            ))

            # Optional tip label
            if label_axes:
                fig.add_trace(go.Scatter3d(
                    x=[end[0]],
                    y=[end[1]],
                    z=[end[2]],
                    mode='text',
                    text=[rf'n_{idx}'],
                    showlegend=False
                ))

        # Origin marker
        fig.add_trace(go.Scatter3d(
            x=[origin[0]],
            y=[origin[1]],
            z=[origin[2]],
            mode='markers',
            marker=dict(size=4, color='white'),
            showlegend=False
        ))

    # ------------------------------------------------------------
    #  Layout
    # ------------------------------------------------------------
    fig.update_layout(
        title=title,
        scene=dict(
            aspectmode='data',
            xaxis=dict(visible=show_axes),
            yaxis=dict(visible=show_axes),
            zaxis=dict(visible=show_axes),
        ),
        margin=dict(l=0, r=0, t=40, b=0),
        legend=dict(itemsizing='constant')
    )
    import plotly.io as pio
    pio.renderers.default = "notebook"

    fig.show()
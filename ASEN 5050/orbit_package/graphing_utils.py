import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
try:
    from orbit_package.ephemeris import Ephemeris
except:
    from ephemeris import Ephemeris
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
import numpy as np
import random
from plotly.subplots import make_subplots
try:
    from orbit_package.orbit_package import Orbit
except:
    from orbit_package import Orbit
orb = Orbit()


import plotly.graph_objects as go
import numpy as np

def get_random_hex_color():
    """Generates a random hexadecimal color code."""
    return f'#{random.randint(0, 0xFFFFFF):06x}'
class Graph:
    def plot_dt(self, dt: list, val: list, name: list) -> None:
        """
        Plot dt vs val with labels using Plotly
        """
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=dt, y=val, mode='lines', name=name[0]
        ))

        fig.update_layout(
            title=f'Graphing {name[0]} Against dt',
            xaxis_title='Time from Initial Condition (s)',
            yaxis_title=f'{name[0]} ({name[1]})',
            showlegend=True
        )
        fig.show()

    def plot_dt_simple(self, dt: list, val: list, name: list) -> None:
        """
        Simplified version of plot_dt
        """
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=dt, y=val, mode='lines', name=name[0]
        ))

        fig.update_layout(
            title=f'Graphing {name[0]} Against dt',
            xaxis_title='Time from Initial Condition (s)',
            yaxis_title=f'{name[0]} ({name[1]})',
            showlegend=True
        )
        fig.show()

    def plot_eph_pos_vector(self, eph, linestyle: str = 'line', gradient: bool = False, earth: bool = False,r_earth:float=6378) -> None:
        """
        Plot in 3D the position vectors for an Ephemeris object using Plotly.
        """
        r_list = eph.all_r()
        t_list = eph.all_dt()

        x_list = [r[0] for r in r_list]
        y_list = [r[1] for r in r_list]
        z_list = [r[2] for r in r_list]

        fig = go.Figure()

        if linestyle.lower() == 'line':
            if gradient:
                fig.add_trace(go.Scatter3d(
                    x=x_list, y=y_list, z=z_list,
                    mode='lines',
                    line=dict(width=4, color=t_list, colorscale='Viridis'),
                    name='Orbit Curve'
                ))
            else:
                fig.add_trace(go.Scatter3d(
                    x=x_list, y=y_list, z=z_list,
                    mode='lines',
                    line=dict(width=4, color='blue'),
                    name='Orbit Curve'
                ))

        elif linestyle.lower() == 'scatter':
            fig.add_trace(go.Scatter3d(
                x=x_list, y=y_list, z=z_list,
                mode='markers',
                marker=dict(size=4, color='red'),
                name='Orbit Measurements'
            ))

        if earth:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            r = r_earth

            xe = r * np.outer(np.cos(u), np.sin(v))
            ye = r * np.outer(np.sin(u), np.sin(v))
            ze = r * np.outer(np.ones(np.size(u)), np.cos(v))

            fig.add_trace(go.Surface(
                x=xe, y=ye, z=ze,
                colorscale=[[0, 'blue'], [1, 'blue']],
                opacity=0.7,
                showscale=False,
                name="Earth"
            ))

            theta = np.linspace(0, 2*np.pi, 500)
            xe_eq = r * np.cos(theta)
            ye_eq = r * np.sin(theta)
            ze_eq = np.zeros_like(theta)

            fig.add_trace(go.Scatter3d(
                x=xe_eq, y=ye_eq, z=ze_eq,
                mode='lines',
                line=dict(color='red', width=5),
                name='Equator'
            ))

        fig.update_layout(
            scene=dict(
                xaxis_title='x [km]',
                yaxis_title='y [km]',
                zaxis_title='z [km]',
                aspectmode = 'data'
            ),
            title="Propagated Orbit Trajectory" if linestyle == 'line' else "Propagated Orbit Position",
            legend=dict(x=0.02, y=0.98)
        )

        fig.show()

    def plot_eph_v_norm(self, ephemeris) -> None:
        """
        Plot velocity norm vs time
        """
        t = ephemeris.all_t()
        v_list = ephemeris.all_v()
        v = [np.linalg.norm(x) for x in v_list]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t, y=v, mode='lines', name='Velocity'
        ))

        fig.update_layout(
            title=f'Velocity vs Time ({ephemeris.data[0][-1]})',
            xaxis_title='Time (UTC)',
            yaxis_title='Velocity (km/s)',
            showlegend=True
        )
        fig.show()

    def plot_eph_r_norm(self, ephemeris) -> None:
        """
        Plot radius norm vs time
        """
        t = ephemeris.all_t()
        r_list = ephemeris.all_r()
        r = [np.linalg.norm(x) for x in r_list]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t, y=r, mode='lines', name='Radius'
        ))

        fig.update_layout(
            title=f'Position Norm vs Time ({ephemeris.data[0][-1]})',
            xaxis_title='Time (UTC)',
            yaxis_title='Radius (km)',
            showlegend=True
        )
        fig.show()

    def plot_eph_nu(self, ephemeris) -> None:
        """
        Plot true anomaly vs time
        """
        t = ephemeris.all_t()
        var_list = ephemeris.all_nu()

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t, y=var_list, mode='lines', name='nu'
        ))

        fig.update_layout(
            title=f'True Anomaly vs Time ({ephemeris.data[0][-1]})',
            xaxis_title='Time (UTC)',
            yaxis_title='True Anomaly (rad)',
            showlegend=True
        )
        fig.show()

    def plot_eph_sma(self, ephemeris) -> None:
        """
        Plot semi-major axis vs time
        """
        t = ephemeris.all_t()
        var_list = ephemeris.all_sma()

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t, y=var_list, mode='lines', name='sma'
        ))

        fig.update_layout(
            title=f'Semi Major Axis vs Time ({ephemeris.data[0][-1]})',
            xaxis_title='Time (UTC)',
            yaxis_title='Semi Major Axis (km)',
            showlegend=True
        )
        fig.show()

    def plot_eph_inc(self, ephemeris) -> None:
        """
        Plot inclination vs time
        """
        t = ephemeris.all_t()
        var_list = ephemeris.all_inc()

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t, y=var_list, mode='lines', name='inc'
        ))

        fig.update_layout(
            title=f'Inclination vs Time ({ephemeris.data[0][-1]})',
            xaxis_title='Time (UTC)',
            yaxis_title='Inclination (rad)',
            showlegend=True
        )
        fig.show()

    def plot_eph_rv(self, eph) -> None:
        """
        Plot r-v space for an Ephemeris
        """
        r = eph.all_r_norm()
        v = eph.all_v_norm()

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=r, y=v, mode='lines', name='R-V Curve'
        ))

        fig.update_layout(
            title='Propagated Orbit R-V Space',
            xaxis_title='Radius (km)',
            yaxis_title='Velocity (km/s)',
            showlegend=True
        )
        fig.show()

    def plot_eph_pos_vector_listeph(self, eph:list, linestyle: str = 'line', gradient: bool = False, earth: bool = False,legend_values:list=None) -> None:
        """
        Plot in 3D the position vectors for a list of eph objects
        """
        fig = go.Figure()
        for i in range(len(eph)):
            if legend_values:
                val = legend_values[i]
            r_list = eph[i].all_r()
            t_list = eph[i].all_dt()

            x_list = [r[0] for r in r_list]
            y_list = [r[1] for r in r_list]
            z_list = [r[2] for r in r_list]

            

            if linestyle.lower() == 'line':
                if gradient:
                    fig.add_trace(go.Scatter3d(
                        x=x_list, y=y_list, z=z_list,
                        mode='lines',
                        line=dict(width=4, color=t_list, colorscale='Viridis'),
                        name=f'Orbit Curve {val}'
                    ))
                else:
                    fig.add_trace(go.Scatter3d(
                        x=x_list, y=y_list, z=z_list,
                        mode='lines',
                        line=dict(width=4, color=get_random_hex_color()),
                        name=f'Orbit Curve {val}'
                    ))

            elif linestyle.lower() == 'scatter':
                fig.add_trace(go.Scatter3d(
                    x=x_list, y=y_list, z=z_list,
                    mode='markers',
                    marker=dict(size=4, color=get_random_hex_color()),
                    name=f'Orbit Curve {val}'
                ))

        if earth:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            r = 6378

            xe = r * np.outer(np.cos(u), np.sin(v))
            ye = r * np.outer(np.sin(u), np.sin(v))
            ze = r * np.outer(np.ones(np.size(u)), np.cos(v))

            fig.add_trace(go.Surface(
                x=xe, y=ye, z=ze,
                colorscale=[[0, 'blue'], [1, 'blue']],
                opacity=0.7,
                showscale=False,
                name="Earth"
            ))

            theta = np.linspace(0, 2*np.pi, 500)
            xe_eq = r * np.cos(theta)
            ye_eq = r * np.sin(theta)
            ze_eq = np.zeros_like(theta)

            fig.add_trace(go.Scatter3d(
                x=xe_eq, y=ye_eq, z=ze_eq,
                mode='lines',
                line=dict(color='red', width=5),
                name='Equator'
            ))

        fig.update_layout(
            scene=dict(
                xaxis_title='x [km]',
                yaxis_title='y [km]',
                zaxis_title='z [km]',
                aspectmode = 'data'
            ),
            title="Propagated Orbit Trajectory" if linestyle == 'line' else "Propagated Orbit Position",
            legend=dict(x=0.02, y=0.98)
            )

        fig.show()

    def plot_eph_diff(self,eph:list,sp:bool=False)-> None:
        """ 
        This function takes in a list of two ephemerides and plots the difference between them
        """
        diff = orb.difference_eph(eph)

        if sp:
            fig = make_subplots(rows=7, cols=1,
                        subplot_titles=("Total Distance", "Delta x", "Delta y", "Delta z","Delta xdot", "Delta ydot", "Delta zdot"))
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dr'],name='dr'),row=1,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dx'],name='dx'),row=2,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dy'],name='dy'),row=3,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dz'],name='dz'),row=4,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dxd'],name='dxd'),row=5,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dyd'],name='dyd'),row=6,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dzd'],name='dzd'),row=7,col=1)
            fig.update_yaxes(title_text="r (km)", row=1, col=1)
            fig.update_yaxes(title_text="x (km)", row=2, col=1)
            fig.update_yaxes(title_text="y (km)", row=3, col=1)
            fig.update_yaxes(title_text="z (km)", row=4, col=1)
            fig.update_yaxes(title_text="xdot (km/s)", row=5, col=1)
            fig.update_yaxes(title_text="ydot (km/s)", row=6, col=1)
            fig.update_yaxes(title_text="zdot (km/w)", row=7, col=1)
            fig.update_layout(
                title='Ephemeris Difference Plot',
                xaxis_title='Time (UTC)',
                autosize=False,
                width=800,  # Set the desired width in pixels
                height=800, # Set the desired height in pixels
                showlegend=True
            )
        else:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=diff['t'], y=diff['dr'], mode='lines', name='Difference Between Ephemerides'
            ))
        
            fig.update_layout(
                title='Ephemeris Difference Plot',
                xaxis_title='Time (UTC)',
                yaxis_title='Distance (km)',
                showlegend=True
            )
        fig.show()

    def plot_eph(self,eph:Ephemeris,include_osc:bool=False,deg:str='DEG')->None:
        """ 
        This will create a grid of plots for analysis of ephemeris
        Args:
            eph: Ephemeris object to plot
            include_osc: If including osculating elements on plots
        """

        fig = make_subplots(rows=5, cols=2,
                        subplot_titles=("Radius", "Eccentricity", "Velocity", "Inclination","True Anomaly","Right Ascension of Ascending Node","Semi Major Axis","Argument of Periapsis","Angular Momentum", "Energy"))
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_r_norm(),name='Distance'),row=1,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_v_norm(),name='Velocity'),row=2,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_nu(deg),name='True Anomaly'),row=3,col=1)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_nu_osc(deg),name='Osculating True Anomaly',mode='markers'),row=3,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma(),name='SMA'),row=4,col=1)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma_osc(),mode='markers',name='Osculating SMA'),row=4,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc(),name='Ecc'),row=1,col=2)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc_osc(),mode='markers',name='Osc Ecc'),row=1,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc(deg),name='Inclination'),row=2,col=2)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc_osc(deg),name='Inclination',mode='markers'),row=2,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_raan(deg),name='RAAN'),row=3,col=2)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_raan_osc(deg),name='Osc RAAN',mode='markers'),row=3,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp(deg),name='Argument Periapsis'),row=4,col=2)
        if include_osc:
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp_osc(deg),name='Osc ArgP',mode='markers'),row=4,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_h(),name='Angular Momentum'),row=5,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_energy(),name='Specific Energy'),row=5,col=2)
        fig.update_yaxes(title_text="Radius (km)", row=1, col=1)
        fig.update_yaxes(title_text="Velocity (km/s)", row=2, col=1)
        fig.update_yaxes(title_text=f"True Anomaly ({deg})",row=3,col=1)
        fig.update_yaxes(title_text=f'Semimajor Axis (km)',row=4,col=1)
        fig.update_yaxes(title_text="Eccentricity ()",row=1,col=2)
        fig.update_yaxes(title_text=f"Inlination ({deg})",row=2,col=2)
        fig.update_yaxes(title_text=f"RAAN ({deg})",row=3,col=2)
        fig.update_yaxes(title_text=f"ArgP ({deg})",row=4,col=2)
        fig.update_yaxes(title_text=f"h (km^2/s)",row=5,col=1)
        fig.update_yaxes(title_text=f"Energy (km^2/s^2)",row=5,col=2)
        

        fig.update_layout(
            title={
            'text': "Ephemeris Features",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1400,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )
        
        fig.show()

    def plot_eph_list(self,ephemerides:list[Ephemeris],include_osc:bool=False,deg:str='DEG',legend_updates:list = None)->None:
        """ 
        This will create a grid of plots for analysis of ephemeris
        Args:
            eph: Ephemeris object to plot
            include_osc: If including osculating elements on plots
        """

        fig = make_subplots(rows=4, cols=2,
                        subplot_titles=("Radius", "Eccentricity", "Velocity", "Inclination","True Anomaly","Right Ascension of Ascending Node","Semi Major Axis","Argument of Periapsis"))
        for i in range(len(ephemerides)):
            eph = ephemerides[i]
            if legend_updates:
                update = legend_updates[i]
            else:
                update = ''
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_r_norm(),name=f'Distance ({update})'),row=1,col=1)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_v_norm(),name=f'Velocity({update})'),row=2,col=1)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_nu(deg),name=f'True Anomaly ({update})'),row=3,col=1)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_nu_osc(deg),name=f'Osculating True Anomaly ({update})',mode='markers'),row=3,col=1)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma(),name=f'SMA ({update})'),row=4,col=1)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma_osc(),mode='markers',name=f'Osculating SMA ({update})'),row=4,col=1)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc(),name=f'Ecc ({update})'),row=1,col=2)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc_osc(),mode='markers',name='Osc Ecc'),row=1,col=2)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc(deg),name=f'Inclination ({update})'),row=2,col=2)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc_osc(deg),name=f'Inclination ({update})',mode='markers'),row=2,col=2)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_raan(deg),name=f'RAAN ({update})'),row=3,col=2)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_raan_osc(deg),name=f'Osc RAAN ({update})',mode='markers'),row=3,col=2)
            fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp(deg),name=f'Argument Periapsis ({update})'),row=4,col=2)
            if include_osc:
                fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp_osc(deg),name=f'Osc ArgP ({update})',mode='markers'),row=4,col=2)

        fig.update_yaxes(title_text="Radius (km)", row=1, col=1)
        fig.update_yaxes(title_text="Velocity (km/s)", row=2, col=1)
        fig.update_yaxes(title_text=f"True Anomaly ({deg})",row=3,col=1)
        fig.update_yaxes(title_text=f'Semimajor Axis (km)',row=4,col=1)
        fig.update_yaxes(title_text="Eccentricity ()",row=1,col=2)
        fig.update_yaxes(title_text=f"Inlination ({deg})",row=2,col=2)
        fig.update_yaxes(title_text=f"RAAN ({deg})",row=3,col=2)
        fig.update_yaxes(title_text=f"ArgP ({deg})",row=4,col=2)

        fig.update_layout(
            title={
            'text': "Ephemeris Comparison",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1400,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )
        
        fig.show()
    
    def plot_eph_cart(self,eph:Ephemeris)->None:
        """ 
        Plot x,y,z,xd,yd,zd
        """
        yd = eph.all_yd()
        fig = make_subplots(rows=3, cols=2,
                        subplot_titles=("x", "xdot", "y", "ydot","z","zdot"))
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_x(),name='x'),row=1,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_y(),name='y'),row=2,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_z(),name='z'),row=3,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_xd(),name='xdot'),row=1,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_yd(),name='ydot'),row=2,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_zd(),name='zdot'),row=3,col=2)

        fig.update_yaxes(title_text="x (km)", row=1, col=1)
        fig.update_yaxes(title_text="y (km)", row=2, col=1)
        fig.update_yaxes(title_text=f"z (km)",row=3,col=1)
        fig.update_yaxes(title_text="xdot",row=1,col=2)
        fig.update_yaxes(title_text=f"ydot",row=2,col=2)
        fig.update_yaxes(title_text=f"zdot",row=3,col=2)

        fig.update_layout(
            title={
            'text': "Cartesian Values over Time",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1400,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )
        fig.show()

    def plot_eph_pos_vector_listpoint(self, eph:Ephemeris,points:list, linestyle: str = 'line', gradient: bool = False, earth: bool = False,legend_values:list=None) -> None:
        """
        Plot in 3D the position vectors for a list of eph objects
        """
        fig = go.Figure()
        t_list = eph.all_dt()

        if linestyle.lower() == 'line':
            if gradient:
                fig.add_trace(go.Scatter3d(
                    x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                    mode='lines',
                    line=dict(width=4, color=t_list, colorscale='Viridis'),
                    name=f'Orbit Curve'
                ))
            else:
                fig.add_trace(go.Scatter3d(
                    x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                    mode='lines',
                    line=dict(width=4, color=get_random_hex_color()),
                    name=f'Orbit Curve'
                ))

        elif linestyle.lower() == 'scatter':
            fig.add_trace(go.Scatter3d(
                x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                mode='markers',
                marker=dict(size=4, color=get_random_hex_color()),
                name=f'Orbit Curve'
            ))

        for i in range(len(points)):
            point = points[i]
            if legend_values:
                val = legend_values[i]
            fig.add_trace(go.Scatter3d(
                x=[float(point[0])], y=[float(point[1])], z=[float(point[2])],
                mode='markers',
                marker=dict(size=4, color=get_random_hex_color()),
                name=f'{val}'
            ))

        if earth:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            r = 6378

            xe = r * np.outer(np.cos(u), np.sin(v))
            ye = r * np.outer(np.sin(u), np.sin(v))
            ze = r * np.outer(np.ones(np.size(u)), np.cos(v))

            fig.add_trace(go.Surface(
                x=xe, y=ye, z=ze,
                colorscale=[[0, 'blue'], [1, 'blue']],
                opacity=0.7,
                showscale=False,
                name="Earth"
            ))

            theta = np.linspace(0, 2*np.pi, 500)
            xe_eq = r * np.cos(theta)
            ye_eq = r * np.sin(theta)
            ze_eq = np.zeros_like(theta)

            fig.add_trace(go.Scatter3d(
                x=xe_eq, y=ye_eq, z=ze_eq,
                mode='lines',
                line=dict(color='red', width=5),
                name='Equator'
            ))

        fig.update_layout(
            scene=dict(
                xaxis_title='x [km]',
                yaxis_title='y [km]',
                zaxis_title='z [km]'
            ),
            title="Propagated Orbit Trajectory" if linestyle == 'line' else "Propagated Orbit Position",
            legend=dict(x=0.02, y=0.98)
            )

        fig.show()

    def plot_eph_pos_vector_listpointeph(self, ephemerides:Ephemeris,points:list, linestyle: str = 'line', gradient: bool = False, earth: bool = False,legend_values_eph:list=None,legend_values_p:list= None) -> None:
        """
        Plot in 3D the position vectors for a list of eph objects and points
        """
        fig = go.Figure()
        

        for i in range(len(ephemerides)):
            eph = ephemerides[i]
            t_list = eph.all_t()
            if legend_values_eph:
                val = legend_values_eph[i]
            if linestyle.lower() == 'line':
                if gradient:
                    fig.add_trace(go.Scatter3d(
                        x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                        mode='lines',
                        line=dict(width=4, color=t_list, colorscale='Viridis'),
                        name=f'Orbit Curve {val}'
                    ))
                else:
                    fig.add_trace(go.Scatter3d(
                        x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                        mode='lines',
                        line=dict(width=4, color=get_random_hex_color()),
                        name=f'Orbit Curve {val}'
                    ))

            elif linestyle.lower() == 'scatter':
                fig.add_trace(go.Scatter3d(
                    x=eph.all_x(), y=eph.all_y(), z=eph.all_z(),
                    mode='markers',
                    marker=dict(size=4, color=get_random_hex_color()),
                    name=f'Orbit Curve {val}'
                ))

        for i in range(len(points)):
            point = points[i]
            if legend_values_p:
                val = legend_values_p[i]
            fig.add_trace(go.Scatter3d(
                x=[float(point[0])], y=[float(point[1])], z=[float(point[2])],
                mode='markers',
                marker=dict(size=4, color=get_random_hex_color()),
                name=f'{val}'
            ))

        if earth:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            r = 6378

            xe = r * np.outer(np.cos(u), np.sin(v))
            ye = r * np.outer(np.sin(u), np.sin(v))
            ze = r * np.outer(np.ones(np.size(u)), np.cos(v))

            fig.add_trace(go.Surface(
                x=xe, y=ye, z=ze,
                colorscale=[[0, 'blue'], [1, 'blue']],
                opacity=0.7,
                showscale=False,
                name="Earth"
            ))

            theta = np.linspace(0, 2*np.pi, 500)
            xe_eq = r * np.cos(theta)
            ye_eq = r * np.sin(theta)
            ze_eq = np.zeros_like(theta)

            fig.add_trace(go.Scatter3d(
                x=xe_eq, y=ye_eq, z=ze_eq,
                mode='lines',
                line=dict(color='red', width=5),
                name='Equator'
            ))

        fig.update_layout(
            scene=dict(
                xaxis_title='x [km]',
                yaxis_title='y [km]',
                zaxis_title='z [km]'
            ),
            title="Propagated Orbit Trajectory" if linestyle == 'line' else "Propagated Orbit Position",
            legend=dict(x=0.02, y=0.98)
            )

        fig.show()

    def plot_eph_pos_vector_listpointeph_2d(self, ephemerides:Ephemeris,points:list, linestyle: str = 'line', gradient: bool = False, earth: bool = False,legend_values_eph:list=None,legend_values_p:list= None) -> None:
        """
        Plot in 3D the position vectors for a list of eph objects and points
        """
        fig = go.Figure()
        

        for i in range(len(ephemerides)):
            eph = ephemerides[i]
            t_list = eph.all_t()
            if legend_values_eph:
                val = legend_values_eph[i]
            if linestyle.lower() == 'line':
                if gradient:
                    fig.add_trace(go.Scatter(
                        x=eph.all_x(), y=eph.all_y(),
                        mode='lines',
                        line=dict(width=4, color=t_list, colorscale='Viridis'),
                        name=f'Orbit Curve {val}'
                    ))
                else:
                    fig.add_trace(go.Scatter(
                        x=eph.all_x(), y=eph.all_y(),
                        mode='lines',
                        line=dict(width=4, color=get_random_hex_color()),
                        name=f'Orbit Curve {val}'
                    ))

            elif linestyle.lower() == 'scatter':
                fig.add_trace(go.Scatter(
                    x=eph.all_x(), y=eph.all_y(),
                    mode='markers',
                    marker=dict(size=4, color=get_random_hex_color()),
                    name=f'Orbit Curve {val}'
                ))

        for i in range(len(points)):
            point = points[i]
            if legend_values_p:
                val = legend_values_p[i]
            fig.add_trace(go.Scatter(
                x=[float(point[0])], y=[float(point[1])],
                mode='markers',
                marker=dict(size=4, color=get_random_hex_color()),
                name=f'{val}'
            ))

        if earth:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            r = 6378

            xe = r * np.outer(np.cos(u), np.sin(v))
            ye = r * np.outer(np.sin(u), np.sin(v))
            ze = np.zeros([len(v),1])

            fig.add_trace(go.Surface(
                x=xe, y=ye,
                colorscale=[[0, 'blue'], [1, 'blue']],
                opacity=0.7,
                showscale=False,
                name="Earth"
            ))

            theta = np.linspace(0, 2*np.pi, 500)
            xe_eq = r * np.cos(theta)
            ye_eq = r * np.sin(theta)

            fig.add_trace(go.Scatter(
                x=xe_eq, y=ye_eq,
                mode='lines',
                line=dict(color='red', width=5),
                name='Equator'
            ))

        fig.update_layout(
            
            xaxis_title='x [km]',
            yaxis_title='y [km]',
            title="Propagated Orbit Trajectory" if linestyle == 'line' else "Propagated Orbit Position",
            legend=dict(x=0.02, y=0.98),
            autosize=False,
            width=800,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            )

        fig.show()

    def plot_eph_groundtrack(self,eph:Ephemeris,mu:float)->None:
        """ 
        This function takes in an ephemeris and creates a graph of the ground track for orbit coverage analysis

        Args:
            eph: Ephemeris to plotgroundtrack of, which requires kep frame
            mu: gravitational parameter of planet
        """
        try:
            if eph.frame == 'ECI_TOD':
                print('Provided ephemeris is in an ECI_TOD frame. Converting to ECF_TOD')


            # Finding r_X using Curtis Algorithm 4.5

        except:
            print('Provided ephemeris has no Keplerian elements. Using Oscullating instead')
            
    def plot_eph_osc(self,eph:Ephemeris,mu:float,deg: str = 'DEG')->None:
        """ 
        This function will take in ephemeris and plot the osculating elements

        Args:
            eph: an ephemeris object to plot
            mu: gravitational parameter

        Returns:
            None
        """

        if eph.data[0][5] is None:
            print('Using Keplers EQ to fill Osc elements')
            eph = orb.fill_eph_osculating(eph,mu)
        
        fig = make_subplots(rows=5, cols=2,
                        subplot_titles=("Radius", "Eccentricity", "Velocity", "Inclination","True Anomaly","Right Ascension of Ascending Node","Semi Major Axis","Argument of Periapsis","Angular Momentum", "Energy"))
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_r_norm(),name='Distance'),row=1,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_v_norm(),name='Velocity'),row=2,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_nu_osc(deg),name='Osculating True Anomaly'),row=3,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma_osc(),name='Osculating SMA'),row=4,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc_osc(),name='Osc Ecc'),row=1,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc_osc(deg),name='Osc Inc'),row=2,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_raan_osc(deg),name='Osc RAAN'),row=3,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp_osc(deg),name='Osc ArgP'),row=4,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_h_osc('norm'),name='Angular Momentum'),row=5,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_energy_osc(),name='Specific Energy'),row=5,col=2)
        fig.update_yaxes(title_text="Radius (km)", row=1, col=1)
        fig.update_yaxes(title_text="Velocity (km/s)", row=2, col=1)
        fig.update_yaxes(title_text=f"True Anomaly ({deg})",row=3,col=1)
        fig.update_yaxes(title_text=f'Semimajor Axis (km)',row=4,col=1)
        fig.update_yaxes(title_text="Eccentricity ()",row=1,col=2)
        fig.update_yaxes(title_text=f"Inlination ({deg})",row=2,col=2)
        fig.update_yaxes(title_text=f"RAAN ({deg})",row=3,col=2)
        fig.update_yaxes(title_text=f"ArgP ({deg})",row=4,col=2)
        fig.update_yaxes(title_text=f"h (km^2/s)",row=5,col=1)
        fig.update_yaxes(title_text=f"Energy (km^2/s^2)",row=5,col=2)
        

        fig.update_layout(
            title={
            'text': "Ephemeris Osculating Features",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1600,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )

        
        fig.show()

    def plot_eph_osc_withJ2(self,eph:Ephemeris,averaged:dict,mu:float,deg: str = 'DEG')->None:
        """ 
        This function will take in ephemeris and plot the osculating elements

        Args:
            eph: an ephemeris object to plot
            mu: gravitational parameter

        Returns:
            None
        """

        ibar = averaged['ibar']
        ibar = [i*180/np.pi for i in ibar]
        Mbar = averaged['Mbar']
        Mbar = [M for M in Mbar]
        Omegabar = averaged['Omegabar']
        Omegabar = [M*180/np.pi + 360 for M in Omegabar]
        omegabar = averaged['omegabar']
        omegabar = [M*180/np.pi for M in omegabar]
        if eph.data[0][5] is None:
            print('Using Keplers EQ to fill Osc elements')
            eph = orb.fill_eph_osculating(eph,mu)
        
        fig = make_subplots(rows=4, cols=2,
                        subplot_titles=("Radius","Argument of Periapsis","Semimajor Axis","Right Ascencsion of the Ascending Node","Eccentricity","Mean Anomaly","Inclination","Mean Motion"))
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_r_norm(),name='Distance'),row=1,col=1)
        # sma subplot (2,1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma_osc(),name='Semimajor Axis'),row=2,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=averaged['abar'],name='Average SMA'),row=2,col=1)
        # ecc subplot (3,1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc_osc(),name='Eccencticity'),row=3,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=averaged['ebar'],name='Average Eccencticity'),row=3,col=1)
        #inc subplot (4,1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc_osc('DEG'),name='Inclination'),row=4,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=ibar,name='Average Inclination'),row=4,col=1)
        # argp subplot (1,2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_argp_osc('DEG'),name='Arg of Perigee'),row=1,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=omegabar,name='Average Argp'),row=1,col=2)
        # raan subplot (2,2)
        fig.add_trace(go.Scatter(x=eph.all_t()[1:-1],y=Omegabar[1:-1],name='Average RAAN'),row=2,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t()[1:-1],y=eph.all_raan_osc('DEG')[1:-1],name='RAAN'),row=2,col=2)
        # M subplot (3,2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_M_osc('DEG'),name='Mean Anomaly'),row=3,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=Mbar,name='Average Mean Anomaly'),row=3,col=2)
        # n subplot (4,2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_n_osc(),name='Mean Motion'),row=4,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=averaged['nbar'],name='Average Mean Motion'),row=4,col=2)
        
        fig.update_yaxes(title_text="R (km)", row=1, col=1)
        fig.update_yaxes(title_text="SMA (km)", row=2, col=1)
        fig.update_yaxes(title_text="ECC", row=3, col=1)
        fig.update_yaxes(title_text=f"INC ({deg})", row=4, col=1)
        fig.update_yaxes(title_text=f"ARGP ({deg})", row=1, col=2)
        fig.update_yaxes(title_text=f"RAAN ({deg})", row=2, col=2)
        fig.update_yaxes(title_text=f"M ({deg})", row=3, col=2)
        fig.update_yaxes(title_text=f"n ({deg}/s)", row=4, col=2)

        fig.update_layout(
            title={
            'text': "Ephemeris Osculating Features with Averaged Results",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1600,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )

        
        fig.show()

    def plot_data(self, y:list, x:list = None, xlab:str=None, ylab:str = None, title:str = None,mode:str='lines',Name:str=None) -> None:
        """
        Plot dt vs val with labels using Plotly
        """
        fig = go.Figure()
        if x is not None:
            fig.add_trace(go.Scatter(
                x=x, y=y, mode=mode,name=Name
            ))
        else:
            x = np.linspace(0,len(y),len(y))
            fig.add_trace(go.Scatter(
                x=x, y=y, mode=mode,name=Name
            ))

        fig.update_layout(
            title=title,
            xaxis_title=xlab,
            yaxis_title=ylab,
            showlegend=True
        )
        fig.show()

    def plot_eph_osc_withsrp(self,eph:Ephemeris,averaged:dict,mu:float,deg: str = 'DEG')->None:
        """ 
        This function will take in ephemeris and plot the osculating elements

        Args:
            eph: an ephemeris object to plot
            mu: gravitational parameter

        Returns:
            None
        """
        abar = averaged['abar']
        ebar = averaged['ebar']
        ibar = averaged['ibar']
        sigmabar = averaged['sigmabar']
        omegatildebar = averaged['omegatildebar']

        M = eph.all_M_osc()
        n = eph.all_n_osc()
        t = eph.all_dt()
        sigma = []
        for k in range(len(eph.data)):
            sigma.append( M[k] - n[k]*t[k])
        
        if eph.data[0][5] is None:
            print('Using Keplers EQ to fill Osc elements')
            eph = orb.fill_eph_osculating(eph,mu)

        omegatilde = eph.all_omegatilde_osc()
        omegatilde = [x*180/np.pi for x in omegatilde]
        
        fig = make_subplots(rows=3, cols=2,
                        subplot_titles=("Radius","Eccentricity","Semimajor Axis","Omega Tilda (omega + Omega)", "Inclination","Sigma"))
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_r_norm(),name='Distance'),row=1,col=1)
        # sma subplot (2,1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_sma_osc(),name='Semimajor Axis'),row=2,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=abar,name='Average SMA'),row=2,col=1)
        # ecc subplot (3,1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_inc_osc('inc'),name='Inclination'),row=3,col=1)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=ibar,name='Average Inclination'),row=3,col=1)
        # argp subplot (1,2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=eph.all_ecc_osc(),name='Eccentricity'),row=1,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=ebar,name='Average Ecc'),row=1,col=2)
        # raan subplot (2,2)
        fig.add_trace(go.Scatter(x=eph.all_t()[1:-1],y=omegatilde,name='Omega tilda'),row=2,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t()[1:-1],y=omegatildebar,name='Avg Omega Tilda'),row=2,col=2)
        # M subplot (3,2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=sigma,name='sigma'),row=3,col=2)
        fig.add_trace(go.Scatter(x=eph.all_t(),y=sigmabar,name='Average sigma'),row=3,col=2)
        
        fig.update_yaxes(title_text="R (km)", row=1, col=1)
        fig.update_yaxes(title_text="SMA (km)", row=2, col=1)
        fig.update_yaxes(title_text=f"INC ({deg})", row=3, col=1)
        fig.update_yaxes(title_text=f"ECC ", row=1, col=2)
        fig.update_yaxes(title_text=f"OMEGA TILDE ({deg})", row=2, col=2)
        fig.update_yaxes(title_text=f"SIGMA ({deg})", row=3, col=2)

        fig.update_layout(
            title={
            'text': "Ephemeris Osculating Features with Averaged Results",
            'font': {'size': 24}
            },
            title_subtitle_text=f"Epoch: {eph.epoch}, Frame: {eph.frame}, Propagation: {eph.data[0][6]}",
            title_subtitle_font={'size': 16, 'color': 'gray'},
            xaxis_title='Time (UTC)',
            autosize=False,
            width=1600,  # Set the desired width in pixels
            height=800, # Set the desired height in pixels
            showlegend=True
        )

        
        fig.show()     
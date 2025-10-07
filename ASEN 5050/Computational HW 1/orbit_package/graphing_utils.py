import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from orbit_package.ephemeris import Ephemeris
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
import numpy as np
import random
from plotly.subplots import make_subplots
from orbit_package.orbit_package import Orbit
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

    def plot_eph_pos_vector(self, eph, linestyle: str = 'line', gradient: bool = False, earth: bool = False) -> None:
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
                zaxis_title='z [km]'
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
            fig = make_subplots(rows=4, cols=1,
                        subplot_titles=("Total Distance", "Delta x", "Delta y", "Delta z"))
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dr'],name='dr'),row=1,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dx'],name='dx'),row=2,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dy'],name='dy'),row=3,col=1)
            fig.add_trace(go.Scatter(x=diff['t'],y=diff['dz'],name='dz'),row=4,col=1)

        else:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=diff['t'], y=diff['dr'], mode='lines', name='Difference Between Ephemerides'
            ))

        fig.update_layout(
            title='Ephemeris Difference Plot',
            xaxis_title='Time (UTC)',
            yaxis_title='Displacement (km)',
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

        fig = make_subplots(rows=4, cols=2,
                        subplot_titles=("Radius", "Eccentricity", "Velocity", "Inclination","True Anomaly","Right Ascension of Ascending Node","Semi Major Axis","Argument of Periapsis"))
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
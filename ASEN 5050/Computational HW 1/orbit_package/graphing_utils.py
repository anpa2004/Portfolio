import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from orbit_package.ephemeris import Ephemeris
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
import numpy as np


class Graph():
    def plot_dt(self,dt:list,val:list,name:list)->None:
        """
        This function takes in an array of dt values and some values corresponding and plots the name of the values

        Args:
            dt: list of delta t values
            val: list of desired thing to graph
            name: labels as strings, of the format:  ['True Anomaly', 'rad'] (ex)
        """

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=dt, y=val, mode='lines', name=name[0]))
        #fig.update_layout(title=f'Graphing {name[0]} Against dt', xaxis_title='Time from Initial Condition (s)', yaxis_title=f'{name[0]} ({name[1]})')
        fig.show()

    def plot_dt_simple(self,dt:list,val:list,name:list)->None:
        """
        This function takes in an array of dt values and some values corresponding and plots the name of the values

        Args:
            dt: list of delta t values
            val: list of desired thing to graph
            name: labels as strings, of the format:  ['True Anomaly', 'rad'] (ex)
        """

        
        plt.plot(dt,val,label=name[0])
        plt.title(f'Graphing {name[0]} Against dt')
        plt.xlabel('Time from Initial Condition (s)')
        plt.ylabel(f'{name[0]} ({name[1]})')
        plt.grid()
        #fig.update_layout(title=f'Graphing {name[0]} Against dt', xaxis_title='Time from Initial Condition (s)', yaxis_title=f'{name[0]} ({name[1]})')
        plt.show()

    def plot_eph_pos_vector(self,eph:Ephemeris,linestyle:str = 'line',gradient:bool = False,earth:bool = False)->None:
        """ 
        This function will be useful to plot in 3d the position vectors for an ehemeris object
        """
        r_list = eph.all_r()
        t_list = eph.all_dt()
        if linestyle.lower() == 'line':
            x_list = []
            y_list = []
            z_list = []
            for r in r_list:
                x_list.append(r[0])
                y_list.append(r[1])
                z_list.append(r[2])

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            if gradient:
                ax.plot3D(x_list, y_list, z_list, label='Orbit Curve',c=t_list)
            else:
                ax.plot3D(x_list, y_list, z_list, label='Orbit Curve')
            ax.set_xlabel('x [km]')
            ax.set_ylabel('y [km]')
            ax.set_zlabel('z [km]')
            ax.set_title('Propagated Orbit Trajectory')
            if earth:
                plt.scatter(0,0,0,c='blue', marker='o',label='Earth')
            ax.legend()
            plt.show()

        elif linestyle.lower() == 'scatter':
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            for r in r_list:
                ax.plot3D(r[0], r[1], r[2], label='Orbit Measurements')
            ax.set_xlabel('x [km]')
            ax.set_ylabel('y [km]')
            ax.set_zlabel('z [km]')
            if earth:
                plt.scatter(0,0,0, c='blue', marker='o',label='Earth')
            ax.set_title('Propagated Orbit Position')
            ax.legend()
            plt.show()
                
    def plot_eph_v_norm(self,ephemeris:Ephemeris)->None:
        """
        This function  plots the norm of velocity vs time
        """        
        t = ephemeris.all_t()
        v_list = ephemeris.all_v()
        v = [np.linalg.norm(x) for x in v_list]
        plt.figure(figsize=(10, 6))
        plt.plot(t,v,label='Velocity')
        plt.grid()
        plt.xlabel('Time (UTC)')
        plt.ylabel('Velocity (km/s)')
        plt.title(f'Velocity vs Time ({ephemeris.data[0][-1]})')
        plt.show()

    def plot_eph_r_norm(self,ephemeris:Ephemeris)->None:
        """
        This function  plots the norm of velocity vs time
        """        
        t = ephemeris.all_t()
        r_list = ephemeris.all_r()
        r = [np.linalg.norm(x) for x in r_list]
        plt.figure(figsize=(10, 6))
        plt.plot(t,r,label='Radius')
        plt.grid()
        plt.xlabel('Time (UTC)')
        plt.ylabel('Radius (km)')
        plt.title(f'Position Norm vs Time ({ephemeris.data[0][-1]})')
        plt.show()

    def plot_eph_nu(self,ephemeris:Ephemeris)->None:
        """
        This function  plots the norm of velocity vs time
        """        
        t = ephemeris.all_t()
        var_list = ephemeris.all_nu()
        plt.figure(figsize=(10, 6))
        plt.plot(t,var_list,label='nu')
        plt.grid()
        plt.xlabel('Time (UTC)')
        plt.ylabel('True Anomaly (Rad)')
        plt.title(f'True Anomaly vs Time ({ephemeris.data[0][-1]})')
        plt.show()

    def plot_eph_sma(self,ephemeris:Ephemeris)->None:
        """
        This function  plots the norm of velocity vs time
        """        
        t = ephemeris.all_t()
        var_list = ephemeris.all_sma()
        plt.figure(figsize=(10, 6))
        plt.plot(t,var_list,label='sma')
        plt.grid()
        plt.xlabel('Time (UTC)')
        plt.ylabel('Semi Major Axis (km)')
        plt.title(f'Semi Major Axis vs Time ({ephemeris.data[0][-1]})')
        plt.show()
    
    def plot_eph_inc(self,ephemeris:Ephemeris)->None:
        """
        This function  plots the norm of velocity vs time
        """        
        t = ephemeris.all_t()
        var_list = ephemeris.all_inc()
        plt.figure(figsize=(10, 6))
        plt.plot(t,var_list,label='inc')
        plt.grid()
        plt.xlabel('Time (UTC)')
        plt.ylabel('Inclination (rad)')
        plt.title(f'Inclination vs Time ({ephemeris.data[0][-1]})')
        plt.show()

    def plot_eph_rv(self,eph:Ephemeris)->None:
        """ 
        This function takes in an ephemeris object and creates an r-v graph
        """
        r = eph.all_r_norm()
        v = eph.all_v_norm()

        plt.plot(r,v)
        plt.xlabel('Radius (km)')  
        plt.ylabel('Velocity (km/s)')  
        plt.title('Propagated orbit R-V Space')
        plt.grid()
        plt.show()

        

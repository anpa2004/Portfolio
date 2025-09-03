import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

class Graph():
    def plot_dt(self,dt:list,val:list,name:list):
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

    def plot_dt_simple(self,dt:list,val:list,name:list):
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

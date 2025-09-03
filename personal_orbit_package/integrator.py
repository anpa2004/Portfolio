""" 
This file will contain the functions for numerical integration to mimic ode45 from MATLAB
"""

import numpy as np
from scipy.integrate import solve_ivp

# An integration function to mimic ode45
def ode45(f,t_span,x0,tval=None,method = 'RK45'):
    if tval is None:
        solution = solve_ivp(f,t_span,x0,method = method)
    else:
        solution = solve_ivp(f,t_span,x0,method=method,t_eval = tval)

    return solution


# Different functions to plug into ode45
mu = 3.986e14 # m^3/s^2

def xd_simple_astro(t,x):
    """ 
    This function calculates xdot for a simple system where an object is orbiting earth and it is assumed
    the mass of the earth is much larger than the object. 

    Args: 
        t: time (not used in xdot calc)
        x: state vector at some time [x,y,z,xd,yd,zd]

    Returns:
        [xd,yd,zd,xdd,ydd,zdd] 
    """

    # Deriving known quantities given current state [x,y,z,xd,yd,zd]
    R = np.array([x[0],x[1],x[2]]) #m
    r = np.linalg.norm(R) #m (scalar distance)
    Rhat = R/r

    # Gravitational attraction to satellite (assuming satellite's is negligible compared to earth) 
    F = -mu/r**2 * Rhat # N
    return np.array([x[3],x[4],x[5],F[0],F[1],F[2]])
import numpy as np
from scipy.interpolate import CubicSpline

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def cot(x):
    return 1/np.tan(x)

def cubic_spline_interpolation(x,y,xq):
    """ 
    This function does interpolation for one or several query points xq

    Args:
        x: list of x data
        y: list of y data
        xq: point or list of points to interpolate at
    
    Returns:
        yq: single value or list of values evaluated at the query point using cubic spline interpolation
    """

    cs = CubicSpline(x,y,bc_type='natural')
    yq = cs(xq)
    return yq

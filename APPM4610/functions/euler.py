import numpy as np

def forward_euler(f,t0,y0,tf,dt):
    """  
    Use a basic Euler Forward step to solve linear first order ODEs
    """
    t_vec = np.arange(t0,tf,dt)

    y = [y0]
    for i in range(len(t_vec)-1):
        yi1 = y[i] + dt*f(t_vec[i+1],y[i])
        y.append(yi1)
    return t_vec,y
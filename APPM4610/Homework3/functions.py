import numpy as np
from scipy.optimize import root_scalar


def forward_euler(f,t0,y0,tf,dt):
    """  
    Use a basic Euler Forward step to solve linear first order ODEs
    imput f = f(t,y)
    """
    N = int(round((tf-t0)/dt))
    t_vec = t0+dt*np.arange(N+1)
    
    y = [y0]
    for i in range(len(t_vec)-1):
        yi1 = y[i] + dt*f(t_vec[i],y[i])
        y.append(yi1)
    return t_vec,y

def ivp_solver(f,t0,y0,tf,dt,method='forward_euler',solver = 'newton',fprime=None,tol=1e-8,nmax=100,neval = False):
    """  
    A general solution function where the user can select a technique to use to solve an IVP
    Inputs are selected based on the user input. 

    Args:
        f: Function f(t,y) to solve
        t0: initial time value
        y0: initial y value
        tf: final integration time
        dt: step size
        method: integration method. defaults to forward_euler. 
            Chose backward_euer, forward_euler, trapezoidal
        solver: type of root finding scheme to use
    """
    n_evaluations = 0
    if method.lower() == 'forward_euler':
        print('Performing forward euler integration')
        N = int(round((tf-t0)/dt))
        t_vec = t0+dt*np.arange(N+1)
        
        y = [y0]
        for i in range(len(t_vec)-1):
            yi1 = y[i] + dt*f(t_vec[i],y[i])
            n_evaluations += 1
            y.append(yi1)
        if not neval:
            return t_vec,y
        if neval:
            return t_vec, y, n_evaluations

    if method.lower() == 'backward_euler':
        print('Performing backward euler')
        N = int(round((tf-t0)/dt))
        t_vec = t0+dt*np.arange(N+1)

        y_vec = [y0]
        for i in range(len(t_vec)-1):
            # Rootfinding for y_{k+1}
            g = lambda y: y - y_vec[-1] - dt*f(t_vec[i+1],y)

            if fprime and solver.lower() == 'newton':
                gp = lambda y: 1 - dt * fprime(t_vec[i+1],y)
                root_results = root_scalar(g,method=solver,x0=y_vec[-1],fprime=gp,xtol=tol,maxiter=nmax)
                n_evaluations += root_results.iterations
            else:
                root_results = root_scalar(g,method=solver,x0=y_vec[-1],xtol=tol,maxiter=nmax)
                n_evaluations += root_results.iterations

            y_vec.append(root_results.root)
        if not neval:
            return t_vec,y_vec
        if neval:
            return t_vec, y_vec, n_evaluations
    
    if method.lower() == 'trapezoidal':
        print('Performing trapezoidal')
        N = int(round((tf-t0)/dt))
        t_vec = t0+dt*np.arange(N+1)

        y_vec = [y0]
        for i in range(len(t_vec)-1):
            g = lambda y: -y+ y_vec[-1] + dt/2*(f(t_vec[i],y_vec[-1]) + f(t_vec[i+1],y))
            if fprime and solver.lower() == 'newton':
                root_results = root_scalar(g,method=solver,x0=y_vec[-1],fprime=fprime,xtol=tol,maxiter=nmax)
                n_evaluations += 2*root_results.iterations
            else:
                root_results = root_scalar(g,method=solver,x0=y_vec[-1],xtol=tol,maxiter=nmax)
                n_evaluations += 2*root_results.iterations
            y_vec.append(root_results.root)
        if not neval:
            return t_vec,y_vec
        if neval:
            return t_vec, y_vec, n_evaluations
    
    if method.lower() == 'rk4':
        print('Performing RK4')

        N = int(round((tf - t0) / dt))
        h = (tf - t0) / N

        ti = t0
        yi = y0

        t = [ti]
        y = [yi]

        for _ in range(N):
            ti = t[-1]
            yi = y[-1]

            k1 = h * f(ti, yi)
            k2 = h * f(ti + h/2, yi + k1/2)
            k3 = h * f(ti + h/2, yi + k2/2)
            k4 = h * f(ti + h, yi + k3)

            t.append(ti + dt)
            y.append(yi + (k1 + 2*k2 + 2*k3 + k4) / 6)
            n_evaluations += 4

        if not neval:
            return np.array(t), np.array(y)
        if neval:
            return np.array(t), np.array(y), n_evaluations

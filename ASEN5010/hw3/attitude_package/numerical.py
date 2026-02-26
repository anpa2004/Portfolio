import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import datetime

class Integrator():
    """ 
    A class filled with several different numerical integration schemes

    ode45: a method for doing numerical integration with a RK45 scheme (technically can be used for any RK scheme)
    
    """
    def ode45(self,f,x0,tspan,other=None,tval=None,method = 'RK45',dense = False,rtol:float=1e-10):
        """ 
        This function is designed to mimic ode45. 

        Args:
            f: some function to integrate. The input requirements must be satisfied by x0 and t, and other
            x0: the initial condition, matching expected input to function f
            tsan: start and end times
            other: any other inputs required for the specific function, such as mu or CD. 
            tval: if integrating at specific timesteps
            method: which integration method to use
            dense: bool of if to return interpolated results

        Returns: 
            t: vector of time values from integration
            x: state of system at each time (each entry is vector depending on input)
        """

        # Setting up a function based on if there are other inputs
        if other:
            def ode_func(t,x):
                return f(t,x,other)
        else:
            def ode_func(t,x):
                return f(t,x)
        
        if tval is None:
            solution = solve_ivp(ode_func,tspan,x0,method = method,dense_output=dense,rtol=rtol)
        else:
            solution = solve_ivp(ode_func,tspan,x0,method=method,t_eval = tval,dense_output=dense,rtol=rtol)

        # Access the results
        t = solution.t
        x = solution.y

        return t,x


    def forward_euler(self,f,t0,y0,tf,dt):
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

    def ivp_solver(self,f,t0,y0,tf,dt,method='forward_euler',solver = 'newton',fprime=None,tol=1e-8,nmax=100):
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

        if method.lower() == 'forward_euler':
            print('Performing forward euler integration')
            N = int(round((tf-t0)/dt))
            t_vec = t0+dt*np.arange(N+1)
            
            y = [y0]
            for i in range(len(t_vec)-1):
                yi1 = y[i] + dt*f(t_vec[i],y[i])
                y.append(yi1)
            return t_vec,y

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
                else:
                    root_results = root_scalar(g,method=solver,x0=y_vec[-1],xtol=tol,maxiter=nmax)

                y_vec.append(root_results.root)
            return t_vec,y_vec
        
        if method.lower() == 'trapezoidal':
            print('Performing trapezoidal')
            N = int(round((tf-t0)/dt))
            t_vec = t0+dt*np.arange(N+1)

            y_vec = [y0]
            for i in range(len(t_vec)-1):
                g = lambda y: -y+ y_vec[-1] + dt/2*(f(t_vec[i],y_vec[-1]) + f(t_vec[i+1],y))
                if fprime and solver.lower() == 'newton':
                    root_results = root_scalar(g,method=solver,x0=y_vec[-1],fprime=fprime,xtol=tol,maxiter=nmax)
                else:
                    root_results = root_scalar(g,method=solver,x0=y_vec[-1],xtol=tol,maxiter=nmax)
                y_vec.append(root_results.root)
            return t_vec,y_vec
        
        
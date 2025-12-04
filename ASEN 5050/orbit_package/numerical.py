import numpy as np
from scipy.integrate import solve_ivp
try:
    from orbit_package.ephemeris import Ephemeris
except:
    from ephemeris import Ephemeris
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
    
    def eph_from_propagation_results(self,t:list,x:list,prop_type:str = '2BODY')-> Ephemeris:
        """ 
        This function takes in the results from propagation and creates an ephemeris object for it

        Args:
            t: list of time values from integration (will be treated as dt)
            x: list of results, assumed to be 2body - [r1,v1]
            prop_type: propagation type- in case x is containing other information
        """
        if prop_type.upper() == '2BODY_DYNBODIES':
            # Propagation from two_body_dynbodies()
            time = [datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(seconds=t0) for t0 in t]
            frame = "ECI_TOD"
            epoch = datetime.datetime.now(datetime.timezone.utc)
            r_list = [v[6:9] for v in x.T]
            v_list = [v[9:12] for v in x.T]
            propagation_type = "NUMERICAL"
            dt = t
            eph = Ephemeris(time,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)
        if prop_type.upper() == '2BODY_NODYN':
            # Propagation from two_body_dynbodies()
            time = [datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(seconds=float(t0)) for t0 in t]
            frame = "ECI_TOD"
            epoch = datetime.datetime.now(datetime.timezone.utc)
            r_list = [v[0:3] for v in x.T]
            v_list = [v[3:6] for v in x.T]
            propagation_type = "NUMERICAL"
            dt = t
            eph = Ephemeris(time,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)
        if prop_type.upper() == '2BODY_J2':
            # Propagation from two_body_dynbodies()
            time = [datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(seconds=float(t0)) for t0 in t]
            frame = "ECI_TOD"
            epoch = datetime.datetime.now(datetime.timezone.utc)
            r_list = [v[0:3] for v in x.T]
            v_list = [v[3:6] for v in x.T]
            propagation_type = "NUMERICAL"
            dt = t
            eph = Ephemeris(time,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)
        if prop_type.upper() == '2BODY_SRP':
            # Propagation from two_body_dynbodies()
            time = [datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(seconds=float(t0)) for t0 in t]
            frame = "ECI_TOD"
            epoch = datetime.datetime.now(datetime.timezone.utc)
            r_list = [v[0:3] for v in x.T]
            v_list = [v[3:6] for v in x.T]
            propagation_type = "NUMERICAL"
            dt = t
            eph = Ephemeris(time,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)
        if prop_type.upper() == '2BODY_SRPJ2':
            # Propagation from two_body_dynbodies()
            time = [datetime.datetime.now(datetime.timezone.utc)+datetime.timedelta(seconds=float(t0)) for t0 in t]
            frame = "ECI_TOD"
            epoch = datetime.datetime.now(datetime.timezone.utc)
            r_list = [v[0:3] for v in x.T]
            v_list = [v[3:6] for v in x.T]
            propagation_type = "NUMERICAL"
            dt = t
            eph = Ephemeris(time,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)
        return eph

        
class ForcingFunction():
    """ 
    This class is a storage place for different ODEs for calculating the derivative state given a current state
    """
    def twobody_nop_dynbodies(self,t,x,other):
        """
        A function to calculate the derivative state for two bodies in 3d space, both with given velocities, masses and positions
        
        Args:
            t: empty value for integrator
            x: a vector [r1,v1,r2,v2] (will be 12 elements)
            other: a dictionary containing m1,m2

        Returns:
            xdot: Derivative as calculated by system dynamics
        """

        m1 = other['m1']
        m2 = other['m2']
        G = other['G']
        r1 = x[0:3]
        v1 = x[3:6]
        r2 = x[6:9]
        v2 = x[9:12]

        # Defining relative vectors
        r12 = r2 - r1
        r = np.linalg.norm(r12)
        rhat = r12/r

        # Calculating Force
        F12 = G*m1*m2/r**2*rhat
        F21 = -F12
        a1 = F12/m1
        a2 = F21/m2

        # Concatenating into a final vector
        xdot = np.concatenate((v1,a1,v2,a2))
        return xdot
    
    def twobody_nop_nodyn(self,t,x,other):
        """ 
        This is a forcing function to numerically integrate a satellite's motion with the following assumptions:
        - The satellite exerts a negligible force on the body it orbits
        - No perturbations
        - The coordinate frame is referenced from the large body's center of mass
        - Everything is a point mass

        Args:
            t: time (unused)
            x: state vector [r1,v1] (will be 6 elements)
            other: dict containing mu

        Returns:
            xdot: Derivative as calculated by system dynamics
        """
        
        r1 = x[0:3]
        v1 = x[3:7]
        mu = other['mu']

        r = np.linalg.norm(r1)
        a = -mu/r**3 * r1

        xdot = np.concatenate((v1,a))
        return xdot
    
    def twobody_j2_nodyn(self,t:float,x:list,other:dict)->float:
        """ 
        This is a forcing function for the propagation of a satellite's initial state under J2 perturbations
        - The satellite exerts a negligible force on the body it orbits
        - No perturbations
        - The coordinate frame is referenced from the large body's center of mass
        - Distances are point masses but earth 
        - J2 harmonic is enough to model the oblateness of the earth

        Args: 
            t: time (unused)
            x: state vector (x,y,z,xd,yd,zd)
            other: dict containing {mu, J2, R, }

        Returns:
            xdot: derivative as calculated by system dynamics
        """

        # Unpacking other:
        mu = other['mu']
        J2 = other['J2']
        R = other['R']

        r1 = x[0:3]
        v1 = x[3:7]

        r = np.linalg.norm(r1)

        x = r1[0]
        y = r1[1]
        z = r1[2]

        Uj2 = -3/2*mu*R**2*J2/r**7*np.array([[(x**2 + y**2 - 4*z**2)*x],[(x**2+y**2 - 4*z**2)*y],[(3*(x**2+y**2) - 2*z**2)*z]])
        a = -mu/r**3 * r1 + np.transpose(Uj2)

        xdot = np.concatenate((v1,a[0]))
        return xdot
    
    def twobody_srp_nodyn(self,t:float,x:list,other:dict)->float:
        """ 
        This is a forcing function for the propagation of a satellite's initial state under Solar Radiation Pressure perturbations
        - The satellite exerts a negligible force on the body it orbits
        - No perturbations
        - The coordinate frame is referenced from the large body's center of mass
        - Distances are point masses but earth 
        - Force is constant from the sun acting in the positive x direction

        Args: 
            t: time (unused)
            x: state vector (x,y,z,xd,yd,zd)
            other: dict containing {mu, g }
                g: scaling factor for solar force

        Returns:
            xdot: derivative as calculated by system dynamics
        """

        # Unpacking other:
        mu = other['mu']
        g = other['g']

        r1 = x[0:3]
        v1 = x[3:7]

        r = np.linalg.norm(r1)

        x = r1[0]
        y = r1[1]
        z = r1[2]

        U = np.array([g,0,0])
        a = -mu/r**3 * r1 + U

        xdot = np.concatenate((v1,a))
        return xdot
    
    def twobody_j2srp_nodyn(self,t:float,x:list,other:dict)->float:
        """ 
        This is a forcing function for the propagation of a satellite's initial state under J2 perturbations
        - The satellite exerts a negligible force on the body it orbits
        - No perturbations
        - The coordinate frame is referenced from the large body's center of mass
        - Distances are point masses but earth 
        - J2 harmonic is enough to model the oblateness of the earth
        - Model SRP as a 1d force that is constant in x direction

        Args: 
            t: time (unused)
            x: state vector (x,y,z,xd,yd,zd)
            other: dict containing {mu, J2, R, g}

        Returns:
            xdot: derivative as calculated by system dynamics
        """

        # Unpacking other:
        mu = other['mu']
        J2 = other['J2']
        R = other['R']
        g = other['g']

        r1 = x[0:3]
        v1 = x[3:7]

        r = np.linalg.norm(r1)

        x = r1[0]
        y = r1[1]
        z = r1[2]

        Uj2 = -3/2*mu*R**2*J2/r**7*np.array([[(x**2 + y**2 - 4*z**2)*x],[(x**2+y**2 - 4*z**2)*y],[(3*(x**2+y**2) - 2*z**2)*z]])
        Usrp = np.array([g,0,0])
        a = -mu/r**3 * r1 + np.transpose(Uj2) + Usrp
        xdot = np.concatenate((v1,a))
        return xdot
        



        
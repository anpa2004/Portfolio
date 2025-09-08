import numpy as np
from scipy.optimize import root_scalar
from orbit_package.math_tools import *
from orbit_package.ephemeris import Ephemeris
import datetime

class Orbit():
    def orbital_elements(self,r_vec:list,v_vec:list,mu:float)->dict:
        """ 
        This function takes in position r, velocity v and gravitational parameter mu and reterns keplerian elements

        Args:
            r_vec: a 3 element vector for position in inertial coordinate frame
            v_vec: a 3 element vector for velocity in inertial coordinate frame
            mu: gravitational parameter in units matching r and v

        Returns:
            ele: a dict containing all the orbital elements at the time of the position and velocity vectors
        """

        # Principle directions
        x_vec = np.array([1,0,0])
        y_vec = np.array([0,1,0])
        z_vec = np.array([0,0,1])

        # Scalar quantities
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)

        # Calculating angular momentum
        h_vec = np.cross(r_vec,v_vec)
        h = np.linalg.norm(h_vec)
        h_hat = h_vec/h

        # Finding inclination
        i = np.arccos(np.dot(h_hat,z_vec))

        # Finding RAAN
        N_vec = np.cross(z_vec,h_vec)
        RAAN = np.arccos(np.dot(N_vec/np.linalg.norm(N_vec),x_vec))

        # Quadrant correction
        if N_vec[1] < 0:
            RAAN = 2*np.pi - RAAN

        # Finding eccentricity
        e_vec = np.cross(v_vec,h_vec)/mu - r_vec/r
        e = np.linalg.norm(e_vec)

        # Finding arg periapsis
        omega = np.arccos(np.dot(e_vec,N_vec)/(e*np.linalg.norm(N_vec)))

        # quadrant correction
        if e_vec[2] <0:
            omega = 2*np.pi - omega

        # Finding true anomaly
        nu = np.arccos(np.dot(e_vec,r_vec)/(e*r))

        # Quadrant correction
        # finding radial velocity
        vr = np.dot(v_vec,r_vec/r)

        if vr<0:
            nu = 2*np.pi - nu

        # Finding sma
        sma = h**2/(mu*(1-e**2))

        # Finding period, mean motion
        T = 2*np.pi*np.sqrt(sma**3/mu)
        n = 2*np.pi/T

        if e<1:
            # Eccentric Anomaly
            E = np.arccos((e+np.cos(nu))/(1+e*np.cos(nu)))

            # Eccentric Anomaly Quadrant Correction
            if nu > np.pi:
                E = 2*np.pi - E

            # Mean Anomaly
            M = E - e*np.sin(E)
        else:
            E = None
            M = None

        if np.abs(i) < np.pi:
            orbit_type = 'prograde'
        else: 
            orbit_type = 'retrograde'

        ele = {"inc": i,"raan":RAAN,"ecc":e,"argp":omega,"nu":nu,"T":T,"sma":sma,"E":E,"M":M,"n":n,"orbit_type":orbit_type}
        return ele

    def true_anomaly_propagation(self,nu:float,e:float,n:float,T:float,t:float)->float:
        """ 
        This function takes in an initial true anomaly and a delta t and propagates 
        the orbit using Kepler's Equation to find the new true anomaly after the dt. 

        Args:
            nu: initial true anomaly
            t: delta t length
            e: orbit eccentricity

        Returns: 
            nuf: final true anomaly (at t0+t)
        """

        E = np.arccos((e+np.cos(nu))/(1+e*np.cos(nu)))
        if nu > np.pi:
            E = 2*np.pi - E
        
        Mi = E - e*np.sin(E)

        # Finding number of revs in the time
        k = t/T
        k = int(k) # finding the integer number of revolutions in that delta t time

        # Finding new M
        Mf = Mi + n*t - 2*k*np.pi
        if t<0:
            Mf = -Mi - n*t + 2*k*np.pi
            Mf = 2*np.pi - Mf

        # Root finding to find new E
        f = lambda E: Mf - E + e*np.sin(E)
        fp = lambda E: -1 + e*np.cos(E)

        Ef = root_scalar(f,method = 'newton',fprime = fp,x0=Mf)
        Ef = Ef.root

        nuf = np.arccos((np.cos(Ef) - e)/(1 - e*np.cos(Ef)))

        if Ef > np.pi:
            nuf = 2*np.pi - nuf
        return nuf,Mf,Ef

    def propagated_elements(self,elements:dict,t:float)->dict:
        """ 
        This function takes in some elements and determines the new orbital elements after the delta t value

        Args:
            elements: a dict from the orbital elements function
            t: a delta t value

        Returns:
            elements_f: final elements at the delta t value
        """

        # using previous function to determine change in nu,M,E
        nu,M,E = self.true_anomaly_propagation(elements['nu'],elements['ecc'],elements['n'],elements['T'],t)

        # Updating elements 
        elements_f = elements
        elements_f['nu'] = nu
        elements_f['M'] = M
        elements_f['E'] = E

        return elements_f
    
    def kep_to_cart(self,elements:dict,mu:float)->list:
        """ 
        This function converts kepelerian elements to cartesian vectors for position and velocity

        Args:
            elements: a dict of keplerian elements
            mu: gravitational parameter

        Returns:
            r_vec: position vector
            v_vec: velocity vector
        """

        # Extracting a few elemends from the elements dict
        a = elements['sma']
        e = elements['ecc']
        E = elements['E']
        nu = elements['nu']
        w = elements['argp']
        i = elements['inc']
        O = elements['raan']

        # Finding distance to central body
        r = a*(1-e*np.cos(E)) # The units of a will dictate the units of r_vec,v_Vec

        # Finding position in orbital frame (z axis is h_vec)
        o = r*np.array([np.cos(nu),np.sin(nu),0])
        od = np.sqrt(a*mu)/r*np.array([-np.sin(E),np.sqrt(1-e**2)*np.cos(E),0])
    
        # Rotating the vectors to their proper frame
        r_vec = []
        v_vec = []

        r_vec.append(o[0]*(cos(w)*cos(O) - sin(w)*cos(i)*sin(O)) - o[1]*(sin(w)*cos(O) + cos(w)*cos(i)*sin(O)))
        r_vec.append(o[0]*(cos(w)*sin(O) + sin(w)*cos(i)*cos(O)) + o[1]*(cos(w)*cos(i)*cos(O) - sin(w)*sin(O)))
        r_vec.append(o[0]*(sin(w)*sin(i)) + o[1]*(cos(w)*sin(i)))

        v_vec.append(od[0]*(cos(w)*cos(O) - sin(w)*cos(i)*sin(O)) - od[1]*(sin(w)*cos(O) + cos(w)*cos(i)*sin(O)))
        v_vec.append(od[0]*(cos(w)*sin(O) + sin(w)*cos(i)*cos(O)) + od[1]*(cos(w)*cos(i)*cos(O) - sin(w)*sin(O)))
        v_vec.append(od[0]*(sin(w)*sin(i)) + od[1]*(cos(w)*sin(i)))

        return r_vec,v_vec

    def propagate_ic(self,r_vec:list,v_vec:list,mu:float,t:float)->list:
        """  
        This function takes in an initial condition (r,v) and propagates for an orbit about an object with
        gravitational parameter mu by time t

        Note: This only returns a value at the time that is away from initial condition by time t (a delta in seconds)

        Args:
            r_vec: initial position vector
            v_vec: initial velocity vector
            mu: gravitational parameter
            t: time difference

        Returns:
            r_f: final position vector
            v_f: final velocity vector
            elementsf: final elements
        """

        # Calculating and propagating the elements
        elements0 = self.orbital_elements(r_vec,v_vec,mu)
        elementsf = self.propagated_elements(elements0,t)
        r_f, v_f = self.kep_to_cart(elementsf,mu)
        return r_f,v_f,elementsf

    def create_propagated_ephemeris(self,r_vec:list,v_vec:list,mu:float,dt:float,epoch:datetime.datetime = None,frame: str = None)-> Ephemeris:    
        """ 
        This uses Kepler's equation to create a propagated ephemeris object
        """

        r_list = []
        v_list = []
        ele_list = []

        propagation_type = 'KEPLER_EQ_PROPAGATION'

        for t in dt:
            r,v,ele = self.propagate_ic(r_vec,v_vec,mu,t)
            r_list.append(r)
            v_list.append(v)
            ele_list.append(ele)

        if not epoch:
            epoch = datetime.datetime.now(datetime.timezone.utc)

        if not frame:
            frame = 'GCF_TOD'
        
        time = []
        for t in dt:
            time.append(epoch + datetime.timedelta(seconds=t))

        eph = Ephemeris(time,elements=ele_list,frame=frame,epoch=epoch,r_vec=r_list,v_vec=v_list,propagation_type = propagation_type,dt=dt)

        return(eph)

    def print_elements(self,ele):
        """ 
        This function is a way of quickly reading out elements after calculating them
        """

        print(f'SMA: {ele['sma']}')
        print(f'RAAN: {ele['raan']}')
        print(f'Arg Periapsis: {ele['argp']}')
        print(f'Inclination: {ele['inc']}')
        print(f'Orbit Period: {ele['T']}')
        print(f'Mean Anomaly: {ele['M']}')
        print(f'Eccentric Anomaly: {ele['E']}')
        print(f'Mean Motion: {ele['n']}')
        print(f'True Anomaly: {ele['nu']}')
        print(f'Orbit Type: {ele['orbit_type']}')
    
    def fill_elements(self,mu:float,inc:float,raan:float,argp:float,ecc:float,rp:float,nu:float)->dict:
        """
        This function takes in minimal information about an orbit and creates the whole elements dict from the information given.
        """

        # Finding semi major axis
        sma = rp/(1-ecc)

        # Finding period
        T = 2*np.pi*sma**(3/2)/np.sqrt(mu)

        n = 2*np.pi/T

        if ecc<1:
            # Eccentric Anomaly
            E = np.arccos((ecc+np.cos(nu))/(1+ecc*np.cos(nu)))

            # Eccentric Anomaly Quadrant Correction
            if nu > np.pi:
                E = 2*np.pi - E

            # Mean Anomaly
            M = E - ecc*np.sin(E)
        else:
            E = None
            M = None

        if np.abs(inc) < np.pi:
            orbit_type = 'prograde'
        else: 
            orbit_type = 'retrograde'

        ele = {"inc": inc,"raan":raan,"ecc":ecc,"argp":argp,"nu":nu,"T":T,"sma":sma,"E":E,"M":M,"n":n,"orbit_type":orbit_type}
        return ele

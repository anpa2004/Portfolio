#### README ####
# This script contains all the  constants and information needed to simulate
# the motion of the Mars orbit. 

###########################################
import numpy as np
from dcms import Attitude_class
from misc_functions import unitize
at = Attitude_class()

class Mars:
    """ 
    Mars class of radius and gravitational parameter, etc
    """
    def __init__(self):
        self.mu = 42828.3 # km^3/s^2
        self.R = 3396.19 # km
        self.rot_period = (1*86400) + (37*60) # seconds
        self.r_gmo = 20424.2 # km

class Mars_Orbit():
    """ 
    Low Mars Orbit parameters
    """ 
    def __init__(self,size):
        mars = Mars()

        # LOW MARS ORBIT
        if size.upper() == 'LMO':
            self.R = mars.R + 400 # km
            self.theta_dot = np.sqrt(mars.mu/self.R**3) # rad/s
            self.raan = np.deg2rad(20) # rad
            self.inc = np.deg2rad(30) # rad
            self.theta0 = np.deg2rad(60) # rad
            self.orbit = 'LMO'

        # GEOSTATIONARY MARS ORBIT
        elif size.upper() == 'GMO':
            self.R = 20424.2 # km
            self.theta_dot = np.sqrt(mars.mu/self.R**3) # rad/s
            self.raan = np.deg2rad(0) # rad
            self.inc = np.deg2rad(0) # rad
            self.theta0 = np.deg2rad(250) # rad
            self.orbit = 'GMO'
        else:
            raise ValueError('Unknown orbit size- expecting LMO or GMO')

    def theta(self,t:float)->float:
        """ 
        Take in time and determine the orbital position for the satellite

        Args:
            t: time to evaluate at

        Returns:
            theta1: theta at time t
            num_revs: number of revolutions since t0
        """

        # Change in angular position (rad)
        dtheta = self.theta_dot*t

        # Correcting for multiple revolutions
        num_revs = dtheta // (np.pi*2)
        dtheta = dtheta % (2*np.pi)

        return dtheta + self.theta0, num_revs
    
    def dcm_in(self,t:float)->list[float]:
        """ 
        Determine the DCM I->N at time t

        Args:
            t: time to evaluate at

        Returns:
            c: dcm at time t
        """

        theta,_ = self.theta(t)

        c = at.dcm_313(self.raan,self.inc,theta)
        return c

    def r_vec(self,t:float)->list[float]:
        """ 
        Determine the r vector at any time in inertial frame

        Args:
            t: time to evaluate at
        
        Returns:
            r in inertial frame
        """
        # Finding dcm
        c = self.dcm_in(t)

        return c.T@np.array([self.R,0,0])

    def v_vec(self,t:float)->list[float]:
        """ 
        Determine the v vector at any time in inertial frame

        Args:
            t: time to evaluate at
        
        Returns:
            r in inertial frame
        """
        # Finding dcm
        c = self.dcm_in(t)

        return c.T@np.array([0,self.R*self.theta_dot,0])
    
    def full_orbital_state(self,t:float,terminal_out:bool=True)->dict:
        """ 
        Return full orbital state for any time

        Args:
            t: time to evaluate at
            terminal_out: boolean of if to print to terminal

        Returns:
            state: dictionary containing raan, inc, theta, r,v,dcm_in
        """
        np.printoptions(precision=4)
        theta,num_revs = self.theta(t)
        c = self.dcm_in(t)
        r = self.r_vec(t)
        v = self.v_vec(t)

        if terminal_out:
            print(f'{self.orbit} at t = {t} seconds:')
            print(f'    r_vec = {r} km')
            print(f'    v_vec = {v} km/s')
            print(f'    number of revolutions since t0: {num_revs}')
            print(f'    theta = {theta} rad = {np.rad2deg(theta)} deg')
            print(f'    raan = {self.raan} rad = {np.rad2deg(self.raan)} deg')
            print(f'    inc = {self.inc} rad = {np.rad2deg(self.inc)} deg')
            print(f'    [IN] = {c}')
            

        return {'theta': theta,'dcm_in':c,'r_vec':r,'v_vec':v}

    def i_r(self,t:float) -> list[float]:
        """
        Determine the radial unit vetor at some time t
        """

        return unitize(self.r_vec(t))
    
    def i_h(self,t:float) -> list[float]:
        """
        Determine the transerse unit vetor at some time t
        """

        return unitize(np.cross(self.r_vec(t),self.v_vec(t)))

    def i_theta(self,t:float) -> list[float]:
        """
        Determine the theta unit vetor at some time t
        """

        return unitize(np.cross(self.i_h(t),self.i_r(t)))
    
    def hill_frame(self,t:float,terminal_out:bool=False)->list[float]: 
        """ 
        Determine the unit vectors and rotation of the hill frame to the inertial
        """
        i_frame = [self.i_r(t),self.i_theta(t),self.i_h(t)]
        n_frame = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
        HN = np.zeros((3,3))

        for i in range(3):
            for j in range(3):
                HN[i,j] = np.dot(i_frame[i],n_frame[j])

        if terminal_out:
            print(f'[HN] at t = {t}:')
            print(f'    [HN] = {HN}')
            print(f'    ir = {i_frame[0]}')
            print(f'    itheta = {i_frame[1]}')
            print(f'    ih = {i_frame[2]}')

        return {'HN':HN,'i_r':i_frame[0],'i_theta':i_frame[1],'i_h':i_frame[2]}
    

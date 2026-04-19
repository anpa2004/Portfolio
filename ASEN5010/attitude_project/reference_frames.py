import numpy as np
from orbits import Mars_Orbit
from misc_functions import unitize, un_tilde


class Reference():
    def __init__(self,size):
        self.orbit = Mars_Orbit(size)

    def sun_pointing_ref(self):
        """ 
        Return information about the sun pointing reference frame
        """
        return {'RN': np.array([[-1,0,0],[0,0,1],[0,1,0]]), 'r1':np.array([1,0,0]), 'r2':np.array([0,1,0]),'r3': np.array([0,0,1]),"omega":np.array([0,0,0])}
    
    def nadir_pointing_ref(self,t:float,terminal_out:bool=False)-> dict:
        """ 
        Take in a time t and generate the nadir pointing frame for the orbit
        """

        # R Frame Basis Vectors
        r1 = - self.orbit.i_r(t)
        r2 = self.orbit.i_theta(t)
        r3 = np.cross(r1,r2)

        r_frame = [r1,r2,r3]
        n_frame = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
        RN = np.zeros((3,3))

        # Creating the DCM
        for i in range(3):
            for j in range(3):
                RN[i,j] = np.dot(r_frame[i],n_frame[j])

        # Creating the angular velocity 
        omega = np.array([0,0,-self.orbit.theta_dot])


        if terminal_out:
            print(f'[RN] at t = {t}:')
            print(f'    [RN] = {RN}')
            print(f'    r1 = {r_frame[0]}')
            print(f'    r2 = {r_frame[1]}')
            print(f'    r3 = {r_frame[2]}')
            print(f'    omega = {omega}')
        return {'RN':RN,'r1':r_frame[0],'r2':r_frame[1],'r3':r_frame[2],"omega":omega}

    def gmo_pointing_ref(self,t:float,dt:float=1e-6,terminal_out:bool=False)-> dict:
        """ 
        Determine the GMO pointing reference frame at some time t

        """

        def vecs2frame(rlmo,rgmo):
            # difference vector delta r
            n3 = np.array([0,0,1])
            dr = rgmo - rlmo
            r1 = - dr/np.linalg.norm(dr)

            # Defining r2,r3
            r2 = unitize(np.cross(dr,n3))
            r3 = np.cross(r1,r2)
            return r1,r2,r3


        if self.orbit.orbit == 'GMO':
            raise ValueError('Cannot use GMO object for GMO reference frame!')
        
        gmo = Mars_Orbit('gmo')
        

        # Reference frame basis vectors
        rLMO = self.orbit.r_vec(t)
        rGMO = gmo.r_vec(t)

        r1,r2,r3 = vecs2frame(rLMO,rGMO)
        
        # forming RN matrix
        RN = np.array([r1.tolist(),r2.tolist(),r3.tolist()])

        ## Determining omega ##################
        # Finding state at t+dt and t-dt
        rLMO_p = self.orbit.r_vec(t + dt)
        rLMO_m = self.orbit.r_vec(t - dt)

        rGMO_p = gmo.r_vec(t + dt)
        rGMO_m = gmo.r_vec(t - dt)

        # Positive matrix
        rp1,rp2,rp3 = vecs2frame(rLMO_p,rGMO_p)
        RN_plus = np.array([rp1.tolist(),rp2.tolist(),rp3.tolist()])

        # Negative matrix
        rm1,rm2,rm3 = vecs2frame(rLMO_m,rGMO_m)
        RN_minus = np.array([rm1.tolist(),rm2.tolist(),rm3.tolist()])

        # build RN_plus and RN_minus
        RN_dot = (RN_plus - RN_minus) / (2*dt)

        # Calculating omega
        omega_tilde = - RN_dot@RN.T
        omega = un_tilde(omega_tilde)

        if terminal_out:
            print(f'[RN] at t = {t}:')
            print(f'    [RN] = {RN}')
            print(f'    r1 = {r1}')
            print(f'    r2 = {r2}')
            print(f'    r3 = {r3}')
            print(f'    omega = {omega}')
        return {'RN':RN,'r1':r1,'r2':r2,'r3':r3,'omega':omega}

import numpy as np
from reference_frames import Reference
from orbits import Mars_Orbit,Mars
from dcms import Attitude_class
from misc_functions import tilde

class Attitude():
    def __init__(self,size,control_schedule = None):
        self.orbit = Mars_Orbit(size)
        self.dcms = Attitude_class()
        self.control_schedule = control_schedule if control_schedule is not None else []

        # LOW MARS ORBIT
        if size.upper() == 'LMO':
            self.mrp_bn_0 = np.array([0.3,-0.4,0.5])                                       # mrp units
            self.omega_bn_0 = np.array([np.deg2rad(1),np.deg2rad(1.75),np.deg2rad(-2.2)])  # rad/s
            self.I_b = np.diag([10,5,7.5])                                                 # kg m^2


        # GEOSTATIONARY MARS ORBIT
        elif size.upper() == 'GMO':
            print('No values for GMO attitude...')
        else:
            raise ValueError('Unknown orbit size- expecting LMO or GMO')
        

    def tracking_error(self,t:float,x:list[float],ref:str = 'sun',terminal_out:bool=False)->dict:
        """ 
        Determine the LMO attitude tracking error 

        Args:
            t: time to evaluate at
            x: state vector
            omega_bn: angular velocity of B with respect to N in INERTIAL coordinates
            RN: rotation matrix R->N

        Returns: 
            errors = {sigma_br,omega_br}, the error in attitude and angular velocity
        """

        # Finding rate of change of the R frame at time t
        reference = Reference(self.orbit.orbit)
        if ref.lower() == 'sun':
            # Sun-pointing reference frame
            reft = ref
            ref = reference.sun_pointing_ref()
            

        elif ref.lower() == 'nadir':
            # Nadir pointing reference frame
            reft = ref
            ref = reference.nadir_pointing_ref(t)

        elif ref.lower() == 'gmo':
            # GMO pointing reference frame
            reft = ref
            ref = reference.gmo_pointing_ref(t)

        else: 
            raise ValueError(f'Error- {ref} is not a recognized reference frame')
        
        # Unpacking reference properties
        r_omega_rn = ref['omega']
        RN = ref['RN']

        # Unpacking state vector
        sigma_bn = x[0:3]
        b_omega_bn = x[3:6]

        # Finding error in attitude
        BN = self.dcms.mrp_to_dcm(sigma_bn)
        BR = BN@RN.T 
        sigma_br = self.dcms.dcm_to_mrp(BR)    

        # Finding the error in omega
        b_omega_rn = BR@r_omega_rn   

        b_omega_br = b_omega_bn - b_omega_rn

        if terminal_out:
            print(f'{reft} pointing frame at time t = {t}:')
            print(f'    BR = {BR}')
            print(f'    sigma_br = {sigma_br}')
            print(f'    b_omega_br = {b_omega_br}')

        return {'BR':BR,'sigma_br':sigma_br,'b_omega_br':b_omega_br}
    
    def rk4_f(self,X:list,t:float,u:list):
        """ 
        Calcualte state derivative given current state and controls
        """
        sigma = X[0:3]
        omega = X[3:6]

        sigma_dot = self.dcms.mrp_der(sigma,omega)
        omega_dot = np.linalg.inv(self.I_b)@(-tilde(omega)@self.I_b@omega + u)

        X_dot = np.hstack((sigma_dot,omega_dot))
        return X_dot

    def scheduled_controller(self,t,X):
        """
        Parse the control_schedule and return a control vector (controlled or not)

        Args:
            t: current time
            X: Current State
        
        Returns
            u: Control vector
        """
        for event in self.control_schedule:
            if t < event['t_end']:
                if t >= event['t_start']:
                    return self.u_pd(t, event['ref'], X, event['K'], event['P'])
                break
        return np.zeros(3)

    def attitude_propagation(self,X0:list[float],tspan:list[float],dt:float=1)->dict[float]:
        """ 
        Propagate the attitude forward given an initial condition and time parameters
        tracking the angular velocity using known system dynamics

        Args:
            X0: Initial State [sigma,omega]
            tspan: Integration Bounds
            u: Control function (function of t and X - u=u(t,X)
            dt: timesetep
        """

        # Parsing Initial Conditions
        sigma0_bn = X0[0:3]
        omega0_bn = X0[3:6]
        I = self.I_b


        # Creating time vector
        N = int(round((tspan[1]-tspan[0])/dt))
        t_vec = tspan[0] + dt*np.arange(N+1)


        # Preallocation 
        X = np.zeros((N+1,6),dtype=float)
        sigma = np.zeros((N+1,3),dtype=float)
        omega = np.zeros((N+1,3),dtype=float)
        u_list = np.zeros((N+1,3),dtype=float)
        X[0] = X0
        sigma[0] = sigma0_bn
        omega[0] = omega0_bn

        H = np.zeros((N+1,3),dtype=float)
        T = np.zeros((N+1,),dtype=float)
        H[0] = I @ omega0_bn
        T[0] = 0.5*omega0_bn.T@I@omega0_bn
        

        # Iterating
        for index,t in enumerate(t_vec[0:-1],0):
            
            # Finding control vector as function of time and position
            u = self.scheduled_controller(t,X[index])
            if np.linalg.norm(u) > 1e3:
                print(f"Large control at t={t}: |u|={np.linalg.norm(u)}")
            u_list[index + 1] = u

            # Performing RK4
            u1 = self.scheduled_controller(t, X[index])
            k1 = dt*self.rk4_f(X[index], t, u1)

            u2 = self.scheduled_controller(t + dt/2, X[index] + k1/2)
            k2 = dt*self.rk4_f(X[index] + k1/2, t + dt/2, u2)

            u3 = self.scheduled_controller(t + dt/2, X[index] + k2/2)
            k3 = dt*self.rk4_f(X[index] + k2/2, t + dt/2, u3)

            u4 = self.scheduled_controller(t + dt, X[index] + k3)
            k4 = dt*self.rk4_f(X[index] + k3, t + dt, u4)
            X[index + 1] = X[index] + 1/6*(k1 + 2*k2 + 2*k3 + k4)

            # Shadow set MRP
            sigma_n1 = X[index + 1][0:3]
            omega_n1 = X[index + 1][3:6]
            if np.linalg.norm(sigma_n1) >= 1:
                sigma_n1 = self.dcms.shadow_mrp(sigma_n1)
                X[index + 1][0:3] = sigma_n1

            # Storing 
            sigma[index + 1] = sigma_n1
            omega[index + 1] = X[index + 1][3:6]
            H[index + 1] = I @ omega_n1
            T[index + 1] = 0.5*omega_n1.T@I@omega_n1

        # Extra helpful parameters
        H_norm = [np.linalg.norm(h) for h in H]
        BN = [self.dcms.mrp_to_dcm(s) for s in sigma]

        return {
            'X': X,
            'sigma': sigma,
            'omega': omega,
            'time': t_vec,
            'u': u_list,
            'T': T,
            'H': H,
            'H_norm':H_norm,
            'BN':BN
        }
        
    def u_pd(self,t,ref,X,K,P):
        """
        Determine the PD control vector for input reference flag and current state

        Args:
            t: current time
            ref: flag of references, matches tracking_error - sun, nadir, gmo
            X: current state
            K: proportional gain
            P: derivative gain

        Returns:
            u: Feedback control vector
        """

        sigma_bn = X[0:3]
        omega_bn = X[3:6]

        errors = self.tracking_error(t,X,ref)

        u = -K*errors['sigma_br'] - P*errors['b_omega_br']
        return u
    
    def rk_integrator(self,X0,t0,tf,u,dt=1):
        """ 
        RK integrate between two points in time with given control u
        """

        # Parsing Initial Conditions
        sigma0_bn = X0[0:3]
        omega0_bn = X0[3:6]
        I = self.I_b


        # Creating time vector
        N = int(round((tf-t0)/dt))
        t_vec = t0 + dt*np.arange(N+1)

        # Preallocation 
        X = np.zeros((N+1,6),dtype=float)
        sigma = np.zeros((N+1,3),dtype=float)
        omega = np.zeros((N+1,3),dtype=float)
        u_list = np.zeros((N+1,3),dtype=float)
        X[0] = X0
        sigma[0] = sigma0_bn
        omega[0] = omega0_bn

        H = np.zeros((N+1,3),dtype=float)
        T = np.zeros((N+1,),dtype=float)
        H[0] = I @ omega0_bn
        T[0] = 0.5*omega0_bn.T@I@omega0_bn

        for index,t in enumerate(t_vec[0:-1],0):
                    
            # Finding control vector as function of time and position
            if u is not None:
                u_vec = u  ### REVISIT WHEN CONTROL IS FURTHER ELABORATED UPON
            else:
                u_vec = np.array([0,0,0])
            u_list[index + 1] = u_vec

            # Performing RK4
            k1 = dt*self.rk4_f(X[index],t,u_vec)
            k2 = dt*self.rk4_f(X[index] + k1/2,t+dt/2,u_vec)
            k3 = dt*self.rk4_f(X[index] + k2/2,t+dt/2,u_vec)
            k4 = dt*self.rk4_f(X[index] + k3,t+dt,u_vec)
            X[index + 1] = X[index] + 1/6*(k1 + 2*k2 + 2*k3 + k4)

            # Shadow set MRP
            sigma_n1 = X[index + 1][0:3]
            omega_n1 = X[index + 1][3:6]
            if np.linalg.norm(sigma_n1) >= 1:
                sigma_n1 = self.dcms.shadow_mrp(sigma_n1)
                X[index + 1][0:3] = sigma_n1

            # Storing 
            sigma[index + 1] = sigma_n1
            omega[index + 1] = X[index + 1][3:6]
            H[index + 1] = I @ omega_n1
            T[index + 1] = 0.5*omega_n1.T@I@omega_n1

        # Extra helpful parameters
        H_norm = [np.linalg.norm(h) for h in H]
        BN = [self.dcms.mrp_to_dcm(s) for s in sigma]
        return {
                    'X': X,
                    'sigma': sigma,
                    'omega': omega,
                    'time': t_vec,
                    'u': u_list,
                    'T': T,
                    'H': H,
                    'H_norm':H_norm,
                    'BN':BN
                }

    
import numpy as np
from dcms import Linear
from scipy.linalg import orth
lin = Linear()

class Dynamics():

    def center_of_gravity_pos_discrete(self,R:list,m:list)->list:
        """  
        Return the center of gravity for a discrete point system
        Finds R_c
        """
        total_mass = sum(m)
        s = 0

        for i in range(len(R)):
            s += 1/total_mass*R[i]*m[i]
        return s
    
    def center_of_gravity_vel_discrete(self,Rdot:list,m:list)->list:
        """  
        Return the center of gravity for a discrete point system
        Finds Rdot_c
        """
        total_mass = sum(m)
        s = 0

        for i in range(len(Rdot)):
            s += 1/total_mass*Rdot[i]*m[i]
        return s
    
    def distance_to_cg_discrete(self,R:list,m:list)->list:
        """ 
        Find the distance to the cg for each point in system
        Finds r
        """
        Rc = self.center_of_gravity_pos_discrete(R,m)
        r = [Ri - Rc for Ri in R]
        return r
    
    def velocity_to_cg_discrete(self,Rd:list,m:list)->list:
        """ 
        Find the velocity relative to cg for each point in system
        Finds rdot
        """
        Rdc = self.center_of_gravity_vel_discrete(Rd,m)
        rd = [Ri - Rdc for Ri in Rd]
        return rd
    
    def kinetic_energy_discrete(self,R:list,Rd:list,m:list)->list:
        """ 
        Return a dict containing rotational kinetic energy, translational and total. 
        """
        # Finding relative position and velocity of each point relative to cg
        r = self.distance_to_cg_discrete(R,m)
        rd = self.velocity_to_cg_discrete(Rd,m)

        # Total mass
        M = sum(m)

        # Translational Kinetic energy
        Rc = self.center_of_gravity_vel_discrete(Rd,m)
        T_trans = 0.5*M*np.dot(Rc,Rc)

        # Rotational Kinetic Energy
        s = 0
        for i in range(len(R)):
            s+= m[i]/2*np.dot(rd[i],rd[i])
        T_rot = s

        return {'total':T_trans + T_rot,'rotational':T_rot, 'Translational':T_trans}

    def angular_momentum_discrete(self,R:list,Rd:list,m:list,rp:list=None,rpd:list=None)->list:
        """ 
        Find the angular momentum of a discrete system about a specified point (rp)
        """

        if rp is None:
            rp = np.zeros(3)
        if rpd is None:
            rpd = np.zeros(3)

        sigma = [r - rp for r in R]
        sigmad = [rd - rpd for rd in Rd]

        Hp = np.zeros(shape=[3])
        for i in range(len(R)):
            Hp += m[i]*np.cross(sigma[i],sigmad[i])

        return Hp
    
    def parallel_axis_theorem(self,Ic:list,Rc:list,M:float)->list:
        """ 
        Compute Io given parameters (Inertia about a point instead of the center of mass)
        """
        return Ic + M*lin.tilde(Rc)@lin.tilde(Rc).T
    
    def similarity_transform(self,IB:list,FB:list)->list:
        """ 
        Perform a change of basis on Ic to change frame (B to F)

        Args:
            ID: Inertia tensor about any point in body frame
            FB: the rotation matrix from B to F

        Returns:
            IB: I represented in F frame
        """

        return FB@IB@FB.T
    
    def principal_inertias(self,I:list)->list:
        """ 
        Calculate and order the principal inertias of an inertia tensor. 

        Args:
            I: Inertia tensor
        
        Returns:
            Ip: Diagonal inertia matrix
        """
        evals,_= np.linalg.eig(I)
        return np.diag(np.flip(np.sort(evals)))
    
    def principle_frame_dcm(self,I:list)->list:
        """ 
        Find the rotation that converts I to the diagonal inertia tensor frame (if it exists for the body)
        This converts to be in descending order of eigenvalues and eigenvectors so I1>I2>I3

        Args:
            I: Inertia matrix

        Return: 
            V: dcm that converts I to diagonal
        """

        # Finding eigenvector matrix
        evals,evecs= np.linalg.eigh(I)
        idx = np.argsort(evals)[::-1]
        V = evecs[:, idx]
        return V.T
import numpy as np
from scipy.optimize import root_scalar


class Linear():
    def tilde(self,v:list)->list:
        """  
        This function calculates the skew matrix representation of a vector v for
        cross product (a x b = tilde{a} * b)

        Args:
            v = vector to create tilde matrix of
        Returns:
            v_tilde = Tilde representation of vector
        """
        try:
            v_tilde = np.array([[0,-v[2], v[1]],
                            [v[2], 0, -v[0]],
                            [-v[1], v[0], 0]])
        except:
            v_tilde = np.array([[0,-v[1][0], v[1][0]],
                            [v[2][0], 0, -v[0][0]],
                            [-v[1][0], v[0][0], 0]])
        return v_tilde
    
    def un_tilde(self,v_tilde:list)->list:
        """   
        Take a tilde matrix and return the corresponding vector
        """
        v = np.array([v_tilde[2,1],v_tilde[0,2],v_tilde[1,0]])
        return v

    def vector(self,x1:float,x2:float,x3:float,shape:str='')->list:
        """  
        This is a simplifying function to create vectors. Args are for a 3x1 vector

        Args: 
            x1: First vector entry
            x2: Second vector entry
            x3: Third vector entry
            shape: 'COL' or 'ROW'
        """
        if shape.upper() == 'COL':
            return np.array([[x1],[x2],[x3]])
        elif shape.upper() == 'ROW':
            return np.array([[x1,x2,x3]])
        else:
            #print('Warning: Unrecognized shape. Expected COL or ROW. Returning (3,)')
            return np.array([x1,x2,x3])

    def vector_to_matrix3(self,v1:list,v2:list,v3:list)->list:
        """ 
        Take in three vectors and formulate a matrix, of the form 
        [V1  V2  V3]
        """
        return np.array([[v1[0],v2[0],v3[0]],[v1[1],v2[1],v3[1]],[v1[2],v2[2],v3[2]]])

    def weighted_outer_product_sum(self,vb:list,vn:list,wk:list):
        """ 
        Do the weighted sum needed to form [K] matrix in Davenport 
        """
        if len(vb) != len(vn) or len(vn)!= len(wk) or len(vb)!=len(wk):
            print('ERROR: Lengths of the vectors or weights are not consistent')

        s = np.zeros(shape=[3,3])
        for i in range(len(vb)):
            s += wk[i]*np.outer(vb[i],vn[i])
        return s
    
    def formulate_k(self,B:list)->list:
        """ 
        Create the k matrix given the sum [B]
        """
        S = B + B.T
        sigma = np.linalg.trace(B)
        Z = np.array([B[1,2] - B[2,1],B[2,0] - B[0,2], B[0,1] - B[1,0]])
        M = S - sigma*np.eye(3)
        return np.array([[sigma,Z[0],Z[1],Z[2]],
                         [Z[0],M[0,0],M[0,1],M[0,2]],
                         [Z[1],M[1,0],M[1,1],M[1,2]],
                         [Z[2],M[2,0],M[2,1],M[2,2]]])

    def sum_diff_measurements(self,vb:list,vn:list)->list:
        """ 
        Find si and di for OLAE
        """
        si = []
        di = []

        for i in range(len(vb)):
            si.append(vb[i] + vn[i])
            di.append(vn[i] - vn[i])
        return si,di
    
lin = Linear()     
class Attitude():
    
    def R1(self,t:float,deg:bool=False)->list:
        """  
        This a rotation about the first axis by angle t 

        Args: 
            t: angle to rotate by
            deg: boolean to use degrees. True = degrees, False = radians (default)
        """
        if deg:
            t = t*np.pi/180
        M = np.array([[1,0,0],[0,np.cos(t),np.sin(t)],[0,-np.sin(t),np.cos(t)]])
        return M

    def R2(self,t:float,deg:bool=False)->list:
        """  
        This a rotation about the second axis by angle t 

        Args: 
            t: angle to rotate by
            deg: boolean to use degrees for user input. True = degrees, False = radians (default)
        """
        if deg:
            t = t*np.pi/180
        M = np.array([[np.cos(t),0,-np.sin(t)],[0,1,0],[np.sin(t),0,np.cos(t)]])
        return M

    def R3(self,t:float,deg:bool=False)->list:
        """  
        This a rotation about the first axis by angle t 

        Args: 
            t: angle to rotate by
            deg: boolean to use degrees. True = degrees, False = radians (default)
        """
        if deg:
            t = t*np.pi/180
        M = np.array([[np.cos(t),np.sin(t),0],[-np.sin(t),np.cos(t),0],[0,0,1]])
        return M

    def dcm_321(self,psi:float,theta:float,phi:float,deg:bool=False)->list:
        """ 
        This creates a DCM for a 3-2-1 Euler angle set
        Args: 
            phi: first angle (not first in succession; angle about first axis. )
            theta: second angle
            psi: third angle
            deg: boolean to use degrees. True = degrees, False = radians (default)
        """
        if deg:
            phi = phi*np.pi/180
            theta = theta*np.pi/180
            psi = psi*np.pi/180

        # Each rotation
        R1 = self.R1(phi)
        R2 = self.R2(theta)
        R3 = self.R3(psi)
        # Composite DCM
        M = R1@R2@R3
        return M
    
    def dcm_313(self,Omega:float,i:float,omega:float,deg:bool=False)->list:
        """ 
        This creates a DCM for a 3-1-3 Euler angle set
        Args: 
            Omega: angle about 1st axis
            i: 1st rotation about 3rd axis
            omega: 2nd rotation about 3rd axis
            deg: boolean to use degrees. True = degrees, False = radians (default)
        """
        if deg:
            Omega = Omega*np.pi/180
            i = i*np.pi/180
            omega = omega*np.pi/180

        c1, s1 = np.cos(Omega), np.sin(Omega)
        c2, s2 = np.cos(i), np.sin(i)
        c3, s3 = np.cos(omega), np.sin(omega)

        C = np.array([
            [ c3*c1 - s3*c2*s1,  c3*s1 + s3*c2*c1,  s3*s2 ],
            [ -s3*c1 - c3*c2*s1, -s3*s1 + c3*c2*c1, c3*s2 ],
            [ s2*s1,            -s2*c1,            c2     ]
        ])
        return C
    
    def dcm_123(self,theta1:float,theta2:float,theta3:float,deg:bool=False)->list:
        """ 
        Create a DCM using a 1-2-3 set of Euler Angles
        """
        return self.R3(theta3,deg=deg)@self.R2(theta2,deg=deg)@self.R1(theta1,deg=deg)

    def inv_321(self,c:list,deg:bool=False):
        """  
        This function will output the list of angles given a DCM for 321 set

        Args: 
            c: dcm to find Euler angles of
            deg: bool of if to return degrees, defaults to false (radians)
        """
        t1 = np.arctan2(c[0,1],c[0,0])
        t2 = -np.arcsin(c[0,2])
        t3 = np.arctan2(c[1,2],c[2,2])

        if deg:
            t1 = t1*180/np.pi
            t2 = t2*180/np.pi
            t3 = t3*180/np.pi
        return t1,t2,t3
    
    def inv_313(self,c:list,deg:bool=False):
        """  
        This function will output the list of angles given a DCM for 313 set

        Args: 
            c: dcm to find Euler angles of
            deg: bool of if to return degrees, defaults to false (radians)
        """
        Omega = np.arctan2(c[2,0],-c[2,1])
        i = np.arccos(c[2,2])
        omega = np.arctan2(c[0,2],c[1,2])

        if deg:
            Omega = Omega*180/np.pi
            i = i*180/np.pi
            omega = omega*180/np.pi
        return Omega, i, omega

    def inv_123(self,c,deg:bool=False):
        """ 
        Invert a 1-2-3 dcm

        """
        # Clamp for numerical safety
        s2 = np.clip(c[2, 0], -1.0, 1.0)
        theta2 = np.arcsin(s2)

        c2 = np.cos(theta2)

        # Check for singularity (gimbal lock)
        if abs(c2) < 1e-8:
            raise ValueError("Singular configuration: cos(theta2) ≈ 0")

        theta1 = np.arctan2(-c[2, 1], c[2, 2])
        theta3 = np.arctan2(-c[1, 0], c[0, 0])
        if deg:
            return theta1*180/np.pi, theta2*180/np.pi, theta3*180/np.pi
        else:
            return theta1, theta2, theta3

    def ea_rate_123(self,theta1,theta2,theta3,omega,deg:bool=False):
        """ 
        Calculate the euler angle rates for a 1-2-3 sequence given a value of omega
        """

        c1,c2,c3 = np.cos(theta1),np.cos(theta2),np.cos(theta3)
        s1,s2,s3 = np.sin(theta1),np.sin(theta2),np.sin(theta3)
        B = 1/c2*np.array([[c3,-s3,0],[c2*s3,c2*c3,0],[-s2*c3,s2*s3,c2]])
        thetad = B@omega
        if deg:
            return thetad[0]*180/np.pi,thetad[1]*180/np.pi,thetad[2]*180/np.pi
        else:
            return thetad[0],thetad[1],thetad[2]


    def crp_to_dcm(self,q:list)->list:
        """  
        Take a standard CRP and create a DCM
        """
        qTq = np.linalg.norm(q)**2
        qqT = np.outer(q,q)
        C = 1/(1+qTq)*((1-qTq)*np.eye(3) + 2*qqT - 2*lin.tilde(q))
        return C
    
    def dcm_to_crp(self, c:list)->list:
        """
        Convert DCM to a CRP
        """ 
        zeta = np.sqrt(np.trace(c) + 1)
        q = 1/(zeta**2)*np.array([c[1,2]-c[2,1],
                                  c[2,0]-c[0,2],
                                  c[0,1]-c[1,0]])
        return q

    def add_crp(self,q1:list,q2:list)->list:
        """  
        Add q1 to q2
        """
        q = (q2 + q1 - np.cross(q2,q1))/(1+np.dot(q2,q1))
        return q
    
    def sub_crp(self,q1:list,q:list)->list:
        """  
        This function finds the relative orientation due to a total and 
        1st q vector. This is effectively subtraction. 
        """
        q2 = (q-q1 + np.cross(q,q1))/(1+np.dot(q,q1))
        return q2

    def crp_der(self,q:list,omega:list)->list:
        """
        Take in a CRP and omega and return the derivative of crp (dot{crp})        
        """
        qdot = 0.5*(np.eye(3)+lin.tilde(q) + np.outer(q,q))@ omega
        return qdot
    
    def prv_from_dcm(self,dcm:list,deg:bool = False)->list:
        """  
        Calculate principle rotation vector from a dcm
        """
        phi = np.arccos(0.5*(np.trace(dcm) - 1))
        ehat = 1/(2*np.sin(phi))*np.array([dcm[1,2] - dcm[2,1],dcm[2,0] - dcm[0,2],dcm[0,1]-dcm[1,0]])
        if deg:
            phi = phi*180/np.pi
        return phi, ehat

    def quat_to_dcm(self,beta:list)->list:
        """   
        Convert a quaternion to a dcm
        """
        b0 = beta[0]
        b1 = beta[1]
        b2 = beta[2]
        b3 = beta[3]

        # Finding elements
        C11 = b0**2 + b1**2 - b2**2 - b3**2
        C12 = 2*(b1*b2 + b0*b3)
        C13 = 2*(b1*b3 - b0*b2)

        C21 = 2*(b1*b2 - b0*b3)
        C22 = b0**2 - b1**2 + b2**2 - b3**2
        C23 = 2*(b2*b3 + b0*b1)

        C31 = 2*(b1*b3 + b0*b2)
        C32 = 2*(b2*b3 - b0*b1)
        C33 = b0**2 - b1**2 - b2**2 + b3**2

        return np.array([[C11, C12, C13],
                        [C21, C22, C23],
                        [C31,C32,C33]])
    
    def dcm_to_quat(self,C:list)->list:
        """  
        Convert a dcm to a quaternion
        """
        b0 = 0.5*np.sqrt(C[0,0] + C[1,1] + C[2,2] + 1)
        b1 = (C[1,2] - C[2,1])/(4*b0)
        b2 = (C[2,0] - C[0,2])/(4*b0)
        b3 = (C[0,1] - C[1,0])/(4*b0)
        return np.array([b0,b1,b2,b3])
    
    def add_quat(self,beta1:list,beta2:list)->list:
        """ 
        Add two quaternions together
        """
        b0 = beta2[0]
        b1 = beta2[1]
        b2 = beta2[2]
        b3 = beta2[3]
        matrix = np.array([[b0,-b1,-b2,-b3],[b1,b0,b3,-b2],[b2,-b3,b0,b1],[b3,b2,-b1,b0]])
        return matrix@beta1
        
    def subtract_quat(self,beta:list,beta1:list)->list:
        """  
        Subtract beta1 from beta
        """
        b0 = beta1[0]
        b1 = beta1[1]
        b2 = beta1[2]
        b3 = beta1[3]
        matrix = np.array([[b0,-b1,-b2,-b3],[b1,b0,b3,-b2],[b2,-b3,b0,b1],[b3,b2,-b1,b0]])
        matrix = np.linalg.inv(matrix)
        return matrix@beta
    
    def shadow_mrp(self,sigma:list)->list:
        """   
        Return the shadow MRP for a given sigma
        """
        s = np.linalg.norm(sigma)
        return -sigma/s**2
    
    def mrp_to_dcm(self,sigma:list)->list:
        """  
        Map an mrp to a dcm
        """
        s = np.linalg.norm(sigma)
        C = np.eye(3) + (8*np.linalg.matrix_power(lin.tilde(sigma),2)-4*(1-s**2)*lin.tilde(sigma))/(1+s**2)**2
        return C
    
    def dcm_to_mrp(self,C:list)->list:
        """
        Return an mrp from a dcm
        """
        zeta = np.sqrt(np.linalg.trace(C) + 1)
        sigma_tilde = (np.transpose(C) - C)/(zeta*(zeta+2))
        
        return lin.un_tilde(sigma_tilde)
    
    def invert_mrp(self,sigma:list)->list:
        """  
        Take an MRP sigma_BN and give sigma_NB
        """
        return -sigma
    
    def add_mrp(self,sigmaFB:list,sigmaBN:list)->list:
        """  
        Add two MRPs together (composite rotation)
        """
        FB = self.mrp_to_dcm(sigmaFB)
        BN = self.mrp_to_dcm(sigmaBN)

        FN = FB@BN
        sigmaFN = self.dcm_to_mrp(FN)
        return sigmaFN
    
    def mrp_der(self,sigma:list,omega:list)->list:
        """  
        Calculate derivative of mrp for a given omega
        """
        s = np.linalg.norm(sigma)
        sigma_dot = 0.25*((1-np.inner(sigma,sigma))*np.eye(3) + 2*lin.tilde(sigma) + 2*np.outer(sigma,sigma))@omega
        return sigma_dot

    def integrate_mrp(self,omega,sigma0:list,tspan:list,dt:float)->list:
        """ 
        Integrate MRPs using shadow sets to avoid singularities
        
        Args:
            omega: Function to calculate omega as a function of time
            sigma0: initial MRP
            tspan: time bounds as a tuple (t0,tf)
            dt: timestep

        """
        # Ensuring initial condition is met
        if np.linalg.norm(sigma0)>1:
            sigma0 = self.shadow_mrp(sigma0)

        # Creating t vector
        N = int(round((tspan[1]-tspan[0])/dt))
        t_vec = tspan[0] + dt*np.arange(N+1)

        t = [tspan[0]]
        sigma = [sigma0]

        # Iterating through the t vals
        while len(t)<=N:
            w = omega(t[-1])
            sigmad = self.mrp_der(sigma[-1],w)
            s = sigma[-1] + dt*sigmad

            if np.linalg.norm(s) > 1:
                s = self.shadow_mrp(s)
            sigma.append(s)
            t.append(t[-1] + dt)
        return t_vec,sigma

    def compare_dcm(self,c1:list,c2:list,deg:bool=False)->list:
        """ 
        Find the principle rotation angle between two DCMs to measure error. Assumed to be given [BN] and [BN]
        """

        BB = c1@c2.T
        phi,_ = self.prv_from_dcm(BB,deg)
        return phi

    def crp_to_quat(self,q:list)->list:
        """ 
        Convert crp to quaternion
        """
        empty_b = np.array([1,q[0],q[1],q[2]])
        beta = 1/np.sqrt(1+np.inner(q,q))*empty_b
        return beta

at = Attitude()
class AttitudeDetermination():
    def triad_method(self,bs:list,bm:list,ns:list,nm:list)->list:
        """ 
        Perform triad method based on two vectors, s and m

        Args:
            bs: reference vector s in body coordinates
            bm: reference vector m in body coordinates
            ns: reference vector s in inertial coordinates
            nm: reference vector m in inertial coordinates
        """

        # Unitizing to machine precision
        bs = bs/np.linalg.norm(bs)
        bm = bm/np.linalg.norm(bm)
        ns = ns/np.linalg.norm(ns)
        nm = nm/np.linalg.norm(nm)

        # Finding t vectors in body coordinates
        bt1 = bs
        bt2 = np.cross(bs,bm)/np.linalg.norm(np.cross(bs,bm))
        bt3 = np.cross(bt1,bt2)/np.linalg.norm(np.cross(bt1,bt2))

        # Finding t vectors in inertial coordinates
        nt1 = ns
        nt2 = np.cross(ns,nm) / np.linalg.norm(np.cross(ns,nm))
        nt3 = np.cross(nt1,nt2) / np.linalg.norm(np.cross(nt1,nt2))

        BT = lin.vector_to_matrix3(bt1,bt2,bt3)
        NT = lin.vector_to_matrix3(nt1,nt2,nt3)

        BN = BT @ np.transpose(NT)
        return BN

    def davenportq_method(self,vb:list,vn:list,wk:list)->list:
        """  
        Take in lists of observations in inertial space and observations 
        in body frame, and a list of weights. 

        Args:
            vb: list of headings in body frame
            vn: list of headings in inertial frame
            wk: list of weights

        Returns:
            dcm: dcm of the orientation
        """

        # Ensuring unit vectors
        vb = [v/np.linalg.norm(v) for v in vb]
        vn = [v/np.linalg.norm(v) for v in vn]
        
        # Finding [K]
        B = lin.weighted_outer_product_sum(vb,vn,wk)
        K = lin.formulate_k(B)

        # Finding Eigenvalues/vectors
        evals,evecs = np.linalg.eig(K)
        evecs = [evecs[:,i] for i in range(len(evecs))]

        # Selecting max eigenvalue and quaternion
        max_index = np.argmax(evals)
        beta = evecs[max_index]

        # Converting to DCM
        dcm = at.quat_to_dcm(beta)
        return dcm
        
    def quest_method(self,vb:list,vn:list,wk:list)->list:
        """ 
        Solve attitude determination problem using QUEST
        """
        # Ensuring unit vectors
        vb = [v/np.linalg.norm(v) for v in vb]
        vn = [v/np.linalg.norm(v) for v in vn]
        
        # Finding [K]
        B = lin.weighted_outer_product_sum(vb,vn,wk)
        K = lin.formulate_k(B)

        # Characteristic equation for rootfinding
        f = lambda s: np.linalg.det(K - s*np.eye(4))
        
        # Rootfinding
        lambda0 = sum(wk)
        root_result = root_scalar(f,x0=lambda0)
        lambda_opt = root_result.root

        # Finding the CRP
        S = B + B.T
        sigma = np.linalg.trace(B)
        Z = np.array([B[1,2] - B[2,1],B[2,0] - B[0,2], B[0,1] - B[1,0]])
        q = np.linalg.inv((lambda_opt + sigma)*np.eye(3) - S)@Z

        # Finding Quaternion
        beta = at.crp_to_quat(q)
        return at.quat_to_dcm(beta)

    def olae_method(self,vb:list,vn:list,wk:list)->list:
        """ 
        Solve Attitude determination problem
        """
        # Ensuring unit vectors
        vb = [v/np.linalg.norm(v) for v in vb]
        vn = [v/np.linalg.norm(v) for v in vn]

        # si and di, d vector
        si,di = lin.sum_diff_measurements(vb,vn)
        d = np.concat(di)

        s_list = []
        w_list = []
        for i in range(len(si)):
            s_list.append(lin.tilde(si[i]))
            w_list.append([wk[i],wk[i],wk[i]])
        w_list = [item for sublist in w_list for item in sublist]

        # Forming S and W
        S = np.vstack(s_list)
        W = np.diag(w_list)

        # Finding q
        q = np.linalg.inv(S.T@W@S)@S.T@W@d
        print(q)
        beta = at.crp_to_quat(q)
        return at.quat_to_dcm(beta)
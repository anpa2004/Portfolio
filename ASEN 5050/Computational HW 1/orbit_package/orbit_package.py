import numpy as np
from scipy.optimize import root_scalar
from orbit_package.math_tools import *
from orbit_package.ephemeris import Ephemeris
import datetime
from orbit_package.math_tools import cubic_spline_interpolation
from numpy.linalg import norm
from scipy.optimize import minimize


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
        if i>0 and i<np.pi:
            N_vec = np.cross(z_vec,h_vec)
            # RAAN = np.arccos(np.dot(N_vec/np.linalg.norm(N_vec),x_vec))
            if np.linalg.norm(N_vec) < 1e-10:
                RAAN = 0
            else:
                RAAN = np.arccos(np.dot(N_vec/np.linalg.norm(N_vec), x_vec))
                if N_vec[1] < 0:
                    RAAN = 2*np.pi - RAAN
        else:
            RAAN = 0
            N_vec = [1,0,0]       

        # Finding eccentricity
        e_vec = np.cross(v_vec,h_vec)/mu - r_vec/r
        e = np.linalg.norm(e_vec)

        # finding radial velocity
        vr = np.dot(v_vec,r_vec/r)

        # Finding arg periapsis
        if e>1e-8:
            omega = np.arccos(np.dot(e_vec,N_vec)/(e*np.linalg.norm(N_vec)))
            # Finding true anomaly
            nu = np.arccos(np.dot(e_vec,r_vec)/(e*r))
            if vr<0:
                nu = 2*np.pi - nu
        else:
            omega = 0
            nu = np.arccos(np.dot(N_vec, r_vec) / (np.linalg.norm(N_vec) * r))
            if r_vec[2] < 0:
                nu = 2*np.pi - nu

        # quadrant correction
        if e_vec[2] <0:
            omega = 2*np.pi - omega        

        # Finding sma
        sma = h**2/(mu*(1-e**2))

        # Finding period, mean motion
        T = 2*np.pi*np.sqrt(sma**3/mu)
        n = 2*np.pi/T

        if e<1:
            # Eccentric Anomaly
            E = 2 * np.arctan(np.sqrt((1-e)/(1+e)) * np.tan(nu/2))

            # Eccentric Anomaly Quadrant Correction
            if nu > np.pi:
                E = 2*np.pi - E

            # Mean Anomaly
            M = E - e*np.sin(E)

            # Hyperbolic anomaly
            H = 0

        else:
            E = None
            H = 2*np.arctanh(np.sqrt((e-1)/(e+1)*np.tan(nu/2)))
            M = -H + e*np.sinh(H)

        # true anomaly equals mean anomaly in circular orbit
        if e == 0:
            nu = M

        if np.abs(i) < np.pi:
            orbit_type = 'prograde'
        else: 
            orbit_type = 'retrograde'
        energy = -mu/r + v**2/2

        ele = {"inc": float(i),"raan":float(RAAN),"ecc":float(e),"argp":float(omega),"nu":float(nu),"T":float(T),"sma":float(sma),"E":float(E),"M":float(M),"H":float(H),"n":float(n),"orbit_type":orbit_type, 'h':h,'h_vec':h_vec,'energy':energy}
        return ele

    def propagated_elements(self, elements: dict, t: float) -> dict:
        """ 
        This function takes in some elements and determines the new orbital elements after the delta t value

        Args:
            elements: a dict from the orbital elements function
            t: a delta t value

        Returns:
            elements_f: final elements at the delta t value
        """
        # copy to avoid mutating caller's dict
        elements_f = elements.copy()

        nu, M, E = self.true_anomaly_propagation(
            elements['nu'], elements['ecc'], elements['n'], elements['T'], t
        )

        # ensure anomalies are normalized to [0,2pi)
        nu = np.mod(nu, 2*np.pi) if np.isfinite(nu) else elements['nu']
        M  = np.mod(M,  2*np.pi) if np.isfinite(M)  else elements['M']
        # E may be None for hyperbolic; keep as-is if not finite

        elements_f['nu'] = nu
        elements_f['M']  = M
        elements_f['E']  = E

        return elements_f

    def true_anomaly_propagation(self, nu: float, e: float, n: float, T: float, t: float):
        """ 
        This function propagates an initial  state using Kepler's Equation. It has been heavily modified to
        make it more robust to numerical instabilities. 

        Args:
            nu: True anomaly
            e: eccentricity
            n: mean motion
            T: Period
            t: change in time

        Returns: 
            nuf: true anomaly after some delta t
            Mf: Mean anomaly after delta t
            Ef: Eccentric anomaly after delta t
        """
        tol = 1e-12

        # Helper: safe clamp
        def safe_arccos(x):
            return np.arccos(np.clip(x, -1.0, 1.0))

        # Elliptic case
        if e < 1.0 - 1e-8:
            # Find E0
            E0 = 2.0 * np.arctan(np.tan(nu/2.0) / np.sqrt((1+e)/(1-e)))
            if not np.isfinite(E0):
                E0 = 0.0

            # Canonicalize to [0, 2pi)
            if E0 < 0:
                E0 += 2*np.pi

            Mi = E0 - e*np.sin(E0)

            # Use modular arithmetic rather than integer k
            Mf = Mi + n*t
            Mf = np.mod(Mf, 2*np.pi)

            # Solve Mf = Ef - e*sin(Ef). Use robust approach.
            f  = lambda Ex: Ex - e*np.sin(Ex) - Mf
            fp = lambda Ex: 1.0 - e*np.cos(Ex)

            # initial guess: use Mf (good for most cases) or E0 for very small t
            x0 = Mf if abs(t) > tol else E0

            Ef = None
            try:
                sol = root_scalar(f, method='newton', fprime=fp, x0=x0, maxiter=60, tol=1e-12)
                Ef = sol.root
                # sometimes Newton returns something outside [0,2pi], normalize
                Ef = np.mod(Ef, 2*np.pi)
                if not np.isfinite(Ef):
                    raise RuntimeError("Non-finite Ef from Newton")
            except Exception:
                # fallback: bracketed bisection in [0, 2pi]
                try:
                    sol = root_scalar(f, method='bisect', bracket=[0.0, 2*np.pi], maxiter=200, xtol=1e-12)
                    Ef = sol.root
                except Exception as ex:
                    # If all fails, return previous anomaly to avoid junk
                    # You may want to log this event
                    # print("Kepler solver failed:", ex)
                    return nu, Mf, E0

            # Convert Ef -> nu using atan2 (continuous)
            sinE = np.sin(Ef)
            cosE = np.cos(Ef)
            denom = 1 - e*cosE
            # guard denom near zero
            if abs(denom) < 1e-16:
                denom = np.copysign(1e-16, denom)
            nuf = np.arctan2(np.sqrt(1 - e**2) * sinE, cosE - e)
            if nuf < 0:
                nuf += 2*np.pi

            # continuity correction: choose the representative of nuf that is closest to input nu
            # so we avoid sudden ±2pi jumps
            delta = nuf - nu
            if delta > np.pi:
                nuf -= 2*np.pi
            elif delta < -np.pi:
                nuf += 2*np.pi

            # Final normalization to [0, 2pi)
            nuf = np.mod(nuf, 2*np.pi)

            return nuf, Mf, Ef

        # Hyperbolic case
        elif e > 1.0 + 1e-8:
            # compute H from nu robustly
            # guard tan(nu/2) large by using atan2 if needed
            tan_half_nu = np.tan(nu/2.0)
            arg = np.sqrt((e-1)/(e+1)) * tan_half_nu
            # If arg >=1 in magnitude arctanh will blow up; handle numerically
            if not np.isfinite(arg):
                arg = np.sign(arg) * (1 - 1e-12)

            H0 = 2.0 * np.arctanh(arg)
            if not np.isfinite(H0):
                H0 = 0.0

            Mi = e*np.sinh(H0) - H0
            Mf = Mi + n*t  # hyperbolic mean-like quantity; no modulo

            f  = lambda Hx: e*np.sinh(Hx) - Hx - Mf
            fp = lambda Hx: e*np.cosh(Hx) - 1

            try:
                sol = root_scalar(f, method='newton', fprime=fp, x0=H0 if abs(t)>tol else Mf, maxiter=100, tol=1e-12)
                Hf = sol.root
                if not np.isfinite(Hf):
                    raise RuntimeError("Non-finite Hf from Newton")
            except Exception:
                # fallback bracket wide
                try:
                    sol = root_scalar(f, method='bisect', bracket=[-1e6, 1e6], maxiter=200)
                    Hf = sol.root
                except Exception:
                    return nu, Mf, H0

            # H -> nu
            nuf = 2.0 * np.arctan2(np.sqrt(e+1)*np.sinh(Hf/2.0), np.sqrt(e-1)*np.cosh(Hf/2.0))
            if nuf < 0:
                nuf += 2*np.pi

            # continuity correction
            delta = nuf - nu
            if delta > np.pi:
                nuf -= 2*np.pi
            elif delta < -np.pi:
                nuf += 2*np.pi
            nuf = np.mod(nuf, 2*np.pi)

            return nuf, Mf, Hf

        # near-parabolic
        else:
            # e extremely close to 1: delegate to universal variable solver (recommended)
            # For now: return input to avoid garbage; better: raise or switch solver
            # raise RuntimeError("e too close to 1. Use universal-variable solver.")
            return nu, elements.get('M', 0.0), None
    
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

    def difference_eph(self,ephemerides:list)->dict:
        """ 
        This function calculates the interpolated differences between two ephemerides
        """
        eph1 = ephemerides[0]
        eph2 = ephemerides[1]

        x1 = eph1.all_x()
        y1 = eph1.all_y()
        z1 = eph1.all_z()
        x2 = eph2.all_x()
        y2 = eph2.all_y()
        z2 = eph2.all_z()
        x1d = eph1.all_x()
        y1d = eph1.all_y()
        z1d = eph1.all_z()
        x2d = eph2.all_x()
        y2d = eph2.all_y()
        z2d = eph2.all_z()
        t1 = eph1.all_t()
        t2 = eph2.all_t()

        # Creating timestamps for interpolation
        t1ts = [t.timestamp() for t in t1]
        t2ts = [t.timestamp() for t in t2]

        # Interpolating r2 to be at t1 timestamps
        x2i = cubic_spline_interpolation(t2ts,x2,t1ts)
        y2i = cubic_spline_interpolation(t2ts,y2,t1ts)
        z2i = cubic_spline_interpolation(t2ts,z2,t1ts)
        x2di = cubic_spline_interpolation(t2ts,x2d,t1ts)
        y2di = cubic_spline_interpolation(t2ts,y2d,t1ts)
        z2di = cubic_spline_interpolation(t2ts,z2d,t1ts)

        dx = [x2i[i]-x1[i] for i in range(len(x1))]
        dy = [y2i[i]-y1[i] for i in range(len(x1))]
        dz = [z2i[i]-z1[i] for i in range(len(x1))]
        dxd = [x2di[i]-x1d[i] for i in range(len(x1))]
        dyd = [y2di[i]-y1d[i] for i in range(len(x1))]
        dzd = [z2di[i]-z1d[i] for i in range(len(x1))]


        dr = []
        for i in range(len(dx)):
            dr.append(np.sqrt(dx[i]**2+dy[i]**2+dz[i]**2))

        diff = {"dx":dx,"dy":dy,"dz":dz,"dr":dr,"t":t1,"dxd":dxd,"dyd":dyd,"dzd":dzd}
        return diff
    
    def fill_eph_osculating(self,eph:Ephemeris,mu:float)->Ephemeris:
        """ 
        This function takes an ephemeris object and calculates the osculating elements and fills them in

        Args: 
            eph: ephemeris object to edit
            mu: gravitational parameter of the orbited body

        Returns:
            eph1: ephemeris object with the osculating elements
        """

        r_list = eph.all_r()
        v_list = eph.all_v()

        elements = []
        for i in range(len(r_list)):
            elements = self.orbital_elements(r_list[i],v_list[i],mu)
            eph.data[i][5] = elements
        
        return eph
    
    def lamberts(self,r1v:list,r2v:list,dt:float,mu:float)->list:
        """
        This function takes in two observed positions of a satellite's orbit and a given 
        transfer time and calculates the velocity vectors at each point to be used in propagation 
        and orbit construction

        Args: 
            r1v: vector of sat at t1
            r2v: vector of sat at t1+t=t2
            dt: transfer time
            mu: gravitational parameter

        Returns:
            v1v: velocity vector at t1
            v2v: velocity vector at t2

        """

        # Finding c
        cv = r2v-r1v
        

        r1 = np.linalg.norm(r1v)
        r2 = np.linalg.norm(r2v)
        c = np.linalg.norm(cv)
        chat = cv/c
        
        # Finding space triangle perimeter
        am = 0.25*(r1+r2+c) # min a
        s = 0.5*(r1+r2+c) # semi perimeter

        # Finding minimum eccentricity
        es = (r2-r1)/c

        # Finding transfer angle
        c_theta = np.dot(r1v,r2v)/(r1*r2)
        theta = np.arccos(c_theta)
        
        # Finding alpha_m
        alpha_m = np.pi
        beta_m = 2*np.arcsin(np.sqrt((s-c)/2))

        # Finding minimum time of flight
        t_m = np.sqrt(mu)**-1*(s**3/8)**0.5*(alpha_m - beta_m + np.sin(beta_m))

        # Finding a for elliptic transfer
        def root_func_ellipse(a,dt=dt,mu=mu):
            """ 
            Function to root solve for the value of a for an elliptic orbit
            Note: a>0
            """
            alpha = 2*np.arcsin(np.sqrt(s/(2*a)))
            beta = 2*np.arcsin(np.sqrt((s-c)/(2*a)))
            f = -dt + np.sqrt(a**3/mu)*(alpha - beta - (sin(alpha) - sin(beta)))
            return f
        
        def root_func_hyperbola(a,dt=dt,mu=mu):
            """ 
            Function to root solve for the value of a for an hyperbola orbit
            Note: a<0
            """
            alpha = 2*np.arcsinh(np.sqrt(s/(-2*a)))
            beta = 2*np.arcsinh(np.sqrt((s-c)/(-2*a)))
            f = -dt + np.sqrt(-a**3/mu)*(np.sinh(alpha)-alpha - np.sinh(beta)+beta)
            return f

        solution = root_scalar(root_func_ellipse,x0=am)
        a_ellipse = solution.root
        solution = root_scalar(root_func_hyperbola,x0=-1)
        a_hyperbola = solution.root

        # Finding alpha and beta
        alpha2_e = np.arcsin(np.sqrt(s/(2*a_ellipse)))
        beta2_e = np.arcsin(np.sqrt((s-c)/(2*a_ellipse)))
        alpha2_h = np.arcsinh(np.sqrt(s/(-2*a_hyperbola)))
        beta2_h = np.arcsinh(np.sqrt((s-c)/(-2*a_hyperbola)))

        # Quadrant corrections for arcsin
        if dt>t_m:
            alpha2_e = 2*np.pi - alpha2_e
        if theta>np.pi and theta<2*np.pi:
            beta_e = -beta_e

        
        # Finding vectors for ellipse
        A_e = np.sqrt(mu/(4*a_ellipse))*cot(alpha2_e)
        B_e = np.sqrt(mu/(4*a_ellipse))*cot(beta2_e)

        v1v_e = (B_e+A_e)*chat + (B_e-A_e)*r1v/r1
        v2v_e = (B_e+A_e)*chat - (B_e-A_e)*r2v/r2

        # Finding vectors for hyperbola
        A_h = np.sqrt(mu/(4*a_hyperbola))*cot(alpha2_h)
        B_h = np.sqrt(mu/(4*a_hyperbola))*cot(beta2_h)

        v1v_h = (B_h+A_h)*chat + (B_h-A_h)*r1v/r1
        v2v_h = (B_h+A_h)*chat - (B_h-A_h)*r2v/r2

        e_ellipse = np.sqrt(1 - (r1*r2*(1-c_theta))/(a_ellipse*(s-c)))
        e_hyperbola = np.sqrt(1 - (r1*r2*(1-c_theta))/(a_hyperbola*(s-c)))

        return v1v_e,v2v_e,v1v_h,v2v_h,a_ellipse,a_hyperbola,e_ellipse,e_hyperbola
    
    def lamberts_tof(self,r1v:list,r2v:list,mu:float)->list:
        """
        This function takes in two observed positions of a satellite's orbit and finds the minimum transfer time
        minimum energy transfer between the two points (am,tm)

        Args: 
            r1v: vector of sat at t1
            r2v: vector of sat at t1+t=t2
            mu: gravitational parameter

        Returns:
            v1v: velocity vector at t1
            v2v: velocity vector at t2
            t_m: time of flight

        """

        # Finding c
        cv = r2v-r1v
        

        r1 = np.linalg.norm(r1v)
        r2 = np.linalg.norm(r2v)
        c = np.linalg.norm(cv)
        chat = cv/c
        
        # Finding space triangle perimeter
        am = 0.25*(r1+r2+c) # min a
        s = 0.5*(r1+r2+c) # semi perimeter

        # Finding transfer angle
        c_theta = np.dot(r1v,r2v)/(r1*r2)
        theta = np.arccos(c_theta)
        
        # Finding alpha_m
        alpha_m = np.pi
        beta_m = 2*np.arcsin(np.sqrt((s-c)/2))

        # Finding minimum time of flight
        t_m = np.sqrt(mu)**-1*(s**3/8)**0.5*(alpha_m - beta_m + np.sin(beta_m))

        # Finding alpha and beta
        alpha2 = np.arcsin(np.sqrt(s/(2*am)))
        beta2 = np.arcsin(np.sqrt((s-c)/(2*am)))
        
        # Finding vectors for ellipse
        A_e = np.sqrt(mu/(4*am))*cot(alpha2)
        B_e = np.sqrt(mu/(4*am))*cot(beta2)

        v1v = (B_e+A_e)*chat + (B_e-A_e)*r1v/r1
        v2v = (B_e+A_e)*chat - (B_e-A_e)*r2v/r2

        return v1v,v2v,t_m

    def lamberts_vallado(self,r0v:list,rv:list,mu:float)->list:
        """ 
        This is Algorithm 56 in Vallado's Astrodynamics textbook, returns the minimum energy transfer properties

        Args:
            r0v: Initial position measurement
            rv: Position measurement some time later
        """
        r0 = norm(r0v)
        r = norm(rv)
        cos_dnu = np.dot(r0v,rv)/(r0*r)
        dnu = np.arccos(cos_dnu)
        c = np.sqrt(r0**2+r**2 - 2*r0*r*cos_dnu)
        s = (r0 + r + c)/2

        #-------------------------------
        # Minimum energy Transfer
        #-------------------------------
        amin = s/2
        pmin = (r*r0)/c * (1 - cos_dnu)
        emin = np.sqrt(1 - 2*pmin/s)

        alphae = np.pi
        betae = 2*np.arcsin(np.sqrt((s - c)/s))

        tmin_aminp = np.sqrt(amin**3/mu)*(alphae + (betae - sin(betae)))
        tmin_aminm = np.sqrt(amin**3/mu)*(alphae - (betae - sin(betae)))

        v0_amin = np.sqrt(mu*pmin)/(r*r0*sin(dnu))*(rv - (1-r/pmin*(1-cos_dnu))*r0v)
        min_energy = {"a": float(amin),"tm":float(tmin_aminm),"tp":float(tmin_aminp),"v0":v0_amin,"e":float(emin),"p":float(pmin)}

        #-------------------------------
        # Minimum eccentricity Transfer
        #-------------------------------
        e_minecc = (r - r0)/c #from P&c
        a_mine = (r*r0*(1-cos_dnu))/((s-c)*(1-e_minecc**2))
        p_mine = a_mine*(1-e_minecc**2)
        tmin_eminp = np.sqrt(a_mine**3/mu)*(alphae + (betae - sin(betae)))
        tmin_eminm = np.sqrt(a_mine**3/mu)*(alphae - (betae - sin(betae)))
        v0_emin = np.sqrt(mu*p_mine)/(r*r0*sin(dnu))*(rv - (1-r/p_mine*(1-cos_dnu))*r0v)
        min_ecc = {"a": float(a_mine),"tm":float(tmin_eminm),"tp":float(tmin_eminp),"v0":v0_emin,"e":float(e_minecc),"p":float(p_mine)}

        #---------------------------------
        # Minimum Time Transfer (parabolic)
        #---------------------------------
        tmin_abs = 1/3*np.sqrt(2/mu)*(s**(3/2) - (s-c)**(3/2))
        ap = None
        ep = 1

        def parabolic_v1(r1_vec, r2_vec, mu, delta_nu, short_way=True):
            r1 = norm(r1_vec)
            r2 = norm(r2_vec)

            # root solve for nu1 from equation (★)
            def f(nu1):
                return r1*(1 + np.cos(nu1)) - r2*(1 + np.cos(nu1 + delta_nu))

            # initial guess for nu1:
            nu1_guess = -delta_nu/2.0
            sol = root_scalar(f, x0=nu1_guess, method='newton')
            nu1 = sol.root

            # compute p (optional check)
            p = r1*(1 + np.cos(nu1))

            # velocity magnitude and components
            vmag = np.sqrt(2*mu/r1)
            gamma = 0.5 * nu1
            vr = vmag * np.sin(gamma)
            vt = vmag * np.cos(gamma)

            # unit vectors (equatorial plane example)
            rhat = r1_vec / r1
            # transverse direction (prograde)
            thetahat = np.array([-rhat[1], rhat[0], 0.0])  # if r vectors are in xy plane
            # if 3D non-equatorial use h-direction or cross-product with k

            v1 = vr * rhat + vt * thetahat
            return v1, p, nu1
        v0p,pp,nu1p = parabolic_v1(r0v,rv,mu,dnu)

        min_time = {"a":ap,"tm":float(tmin_abs),"v0":v0p,"e":float(ep),"p":float(pp)}
        
        return min_energy,min_ecc,min_time
    
    def lamberts_gauss(self,r0v:list,rv:list,dt:float,tm:float,mu:float,e:str,tol:float = 1e-8)->list:
        """ 
        This is outlined in Vallado pg 478, algo 57
        Args:
            r0v: Initial position
            rv: Position some time later
            dt: delta t of transfer
            tm: minimum transfer time

        returns:
            v0: velocity at point r0
            v: velocity at point r
        """
        r = norm(rv)
        r0 = norm(r0v)
        cos_dn = np.dot(r0v,rv)/(r*r0)
        sin_dn = tm*np.sqrt(1-cos_dn**2)
        dn = np.arctan2(sin_dn,cos_dn)

        l = (r0+r)/(4*np.sqrt(r0*r)*cos(dn/2)) - 0.5
        m = mu*dt**2/(2*np.sqrt(r0*r)*cos(dn/2))**3

        y = 1
        y0 = 0
        while y-y0> tol:
            y0 = y
            x1 = m/y0**2 - l
            x2 = 4/3*(1 + 6*x1/5 + 48*x1**2/35 + 680*x1**3/(5*7*9))
            y = 1 + x2*(l+x1)

        if e > 0 and e<1:
            cos_dE2 = 1 - 2*x1
            p = (r*r0*(1 - cos_dn))/(r0 + r - 2*np.sqrt(r0*r)*cos(dn/2)*cos_dE2)

        elif e>1:
            cosh_dH2 = 1 - 2*x1
            p = (r*r0*(1 - cos_dn))/(r0 + r - 2*np.sqrt(r0*r)*cos(dn/2)*cosh_dH2)

        f = 1 - r/p*(1-cos_dn)
        g = r*r0*sin_dn/(np.sqrt(mu*p))
        fd = np.sqrt(1/p)*np.tan(dn/2)*((1-cos_dn)/p - 1/r - 1/r0)
        gd = 1 - r0/p*(1-cos_dn)
        v0 = (rv - f*r0v)/g
        v = (gd*rv - r0)/g

        return v0,v

    def lamberts_all_tofs(self,r1v,r2v,mu):
        """ 
        This is a full collection of Lambert's work. It takes in two position vectos
        and a time of flight. It works out all of the possible orbits between the two 
        points and determines whether the TOF is for an elliptic orbit (short or long)
        a hyperbolic orbit, or a parabolic orbit. IT returns all the velocity and information
        for these orbits to use in creating plots and analysis, but it will return
        mainly the information for the actual orbit. 

        Args:
            r1v: Position Vector 1
            r2v: Position Vector 2
            tof_input: Time of Flight
            mu: Gravitation parameter

        Returns:

        """

        r1 = norm(r1v)
        r2 = norm(r2v)
        c = norm(r2v - r1v)
        theta_short = np.arccos(np.dot(r1v,r2v)/(r1*r2))
        theta_long = 2*np.pi - theta_short

        s = 0.5*(r1 + r2 + c)
        am = 0.5*s

        # # Finding orbit transfer times
        # def root_func_ellipse(a,dt,mu):
        #     """ 
        #     Function to root solve for the value of a for an elliptic orbit
        #     Note: a>0
        #     """
        #     alpha = 2*np.arcsin(np.sqrt(s/(2*a)))
        #     beta = 2*np.arcsin(np.sqrt((s-c)/(2*a)))
        #     f = -dt + np.sqrt(a**3/mu)*(alpha - beta - (sin(alpha) - sin(beta)))
        #     return f
        
        # def root_func_hyperbola(a,dt,mu):
        #     """ 
        #     Function to root solve for the value of a for an hyperbola orbit
        #     Note: a<0
        #     """
        #     alpha = 2*np.arcsinh(np.sqrt(s/(-2*a)))
        #     beta = 2*np.arcsinh(np.sqrt((s-c)/(-2*a)))
        #     f = -dt + np.sqrt(-a**3/mu)*(np.sinh(alpha)-alpha - np.sinh(beta)+beta)
        #     return f

        # solution = root_scalar(root_func_ellipse,x0=am)
        # a_ellipse = solution[0]
        # solution = root_scalar(root_func_hyperbola,x0=-1)
        # a_hyperbola = solution[0]

        a_ellipse = am
        a_hyperbola = -am

        # Elliptic transfer time
        alpha_short = 2*np.arcsin(np.sqrt(s/(2*a_ellipse)))
        alpha_long = 2*np.pi - alpha_short
        beta = 2*np.arcsin((s-c)/(2*a_ellipse))

        t_elliptic_short = np.sqrt(a_ellipse**3/mu)*((alpha_short - beta) - (sin(alpha_short) - sin(beta)))
        t_elliptic_long = np.sqrt(a_ellipse**3/mu)*((alpha_long - beta) - (sin(alpha_long) - sin(beta)))
        
        # Minimum energy transfer time
        alpha_m = np.pi
        beta_m = 2*np.sqrt((s-c)/(2*am))

        t_min_energy = np.sqrt(am**3/mu)*((alpha_m - beta_m) - (sin(alpha_m) - sin(beta_m)))

        # Parabolic time of flight
        D = (s - c)/s
        t_parabolic = 1/3*np.sqrt(2/mu)*((r1+r2+c)**(3/2) - (r1+r2 - c)**(3/2))/np.sqrt(r1*r2)

        # Hyperbolic time of flight
        alpha_hyp = 2*np.arcsinh(np.sqrt(s/(-2*a_hyperbola)))
        beta_hyp = 2*np.arcsinh(np.sqrt((s-c)/(-2*a_hyperbola)))

        t_hyperbolic = np.sqrt(-a_hyperbola**3/mu)*((np.sinh(alpha_hyp) - np.sinh(beta_hyp)) - (alpha_hyp - beta_hyp))

        return {"t_elliptic_long":t_elliptic_long,"t_elliptic_short":t_elliptic_short,"t_min_energy":t_min_energy,"t_parabolic":t_parabolic,"t_hyperbolic":t_hyperbolic}
    
    def lamberts_min_e(self,r1v:list,r2v:list,mu:float)-> list:
        """ 
        This function solves Lambert's to find the time of flight and eccentricity of the minimum eccentricity orbit

        Args:
            r1v: position estimate 1
            r2v: position estimate 2
            mu: gravitational parameter

        Returns: 
            vals: dict of information of min eccentricity ellipse
        
        """

        r1 = norm(r1v)
        r2 = norm(r2v)
        cos_dnu = np.dot(r1v,r2v)/(r1*r2)
        dnu = np.arccos(cos_dnu)
        c = np.sqrt(r1**2+r2**2 - 2*r1*r2*cos_dnu)
        s = (r1 + r2 + c)/2

        amin = s/2
        pmin = (r1*r2)/c * (1 - cos_dnu)

        # Optimization routine to find minimum e and corresponding a
        def min_e(a):
            return np.sqrt(1 - (r1*r2*(1-cos_dnu))/(a*(s-c)))
        
        a0 = amin + 0.01
        result = minimize(min_e,a0,method='BFGS')
        emin = result.fun
        amin = result.x

        alphae = np.pi
        betae = 2*np.arcsin(np.sqrt((s - c)/s))

        tmin_eminp = np.sqrt(amin**3/mu)*(alphae + (betae - sin(betae)))
        tmin_eminm = np.sqrt(amin**3/mu)*(alphae - (betae - sin(betae)))

        v0 = np.sqrt(mu*pmin)/(r1*r2*sin(dnu))*(r2v - (1-r1/pmin*(1-cos_dnu))*r1v)

        return {'amin': amin,"emin":emin,"tmin_eminp":tmin_eminp,"tmin_eminm":tmin_eminm,"v0":v0}


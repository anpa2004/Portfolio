import numpy as np
#from orbit_package import orbit
try:
    from orbit_package.math_tools import *
except:
    from math_tools import *
import datetime

#orb = orbit()
class Ephemeris:
    def __init__(self, t: list, r_vec: list = None, v_vec: list = None,
                 elements: dict = None, frame: str = None,
                 epoch: datetime.datetime = None, propagation_type: str = None,
                 dt: list = None,osculating:list = None,perturbed:list=None):

        if not propagation_type:
            propagation_type = 'NO_PROPAGATION_TYPE'

        data = []
        for i in range(len(t)):
            sub_list = [
                t[i],
                elements[i] if elements is not None else None,
                r_vec[i] if r_vec is not None else None,
                v_vec[i] if v_vec is not None else None,
                dt[i] if dt is not None else None,
                osculating[i] if osculating is not None else None,
                perturbed[i] if perturbed is not None else None,
                propagation_type
            ]
            data.append(sub_list)

        self.data = data
        self.epoch = epoch
        self.frame = frame

    def all_r(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[2])
        return r_list
    
    def all_r_norm(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[2])
        r_list = [np.linalg.norm(r) for r in r_list]
        return r_list
    
    def all_v(self):
        """ 
        This function returns a list of every velocity vector contained in the ephemeris object
        """
        v_list = []
        for dat in self.data:
            v_list.append(dat[3])
        return v_list

    def all_v_norm(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        v_list = []
        for dat in self.data:
            v_list.append(dat[3])
        v_list = [np.linalg.norm(v) for v in v_list]
        return v_list

    def all_ele(self):
        """ 
        This function returns a list of every element with time contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1])
        return ele_list
    
    def all_t(self):
        """ 
        This function returns a list of every time stamo contained in the ephemeris object
        """
        t_list = []
        for dat in self.data:
            t_list.append(dat[0])
        return t_list
    
    def all_dt(self):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        t_list = []
        for dat in self.data:
            t_list.append(dat[4])
        return t_list
    
    def all_nu(self,unit:str='RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['nu'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list
    
    def all_sma(self):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['sma'])
        return ele_list
    
    def all_ecc(self):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['ecc'])
        return ele_list
    
    def all_inc(self,unit:str='RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['inc'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list
    
    def all_raan(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['raan'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list
    
    def all_argp(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[1]['argp'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list
    
    def all_xd(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[3][0])
        return r_list
    
    def all_yd(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[3][1])
        return r_list
    
    def all_zd(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[3][2])
        return r_list
    
    def all_x(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[2][0])
        return r_list
    
    def all_y(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[2][1])
        return r_list
    
    def all_z(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[2][2])
        return r_list
    
    def all_nu_osc(self,unit:str='RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['nu'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list
    
    def all_h(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[1]['h'])
        return r_list
    
    def all_h_osc(self,style:str = 'vec'):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        if style.lower() == 'vec':
            r_list = []
            for dat in self.data:
                r_list.append(dat[5]['h_vec'])
            return r_list
        if style.lower() == 'norm':
            r_list = []
            for dat in self.data:
                r_list.append(dat[5]['h'])
            return r_list
         
    def all_energy(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[1]['energy'])
        return r_list
    
    def all_energy_osc(self):
        """ 
        This function returns a list of every position vector contained in the ephemeris object
        """
        r_list = []
        for dat in self.data:
            r_list.append(dat[5]['energy'])
        return r_list
     
    def all_sma_osc(self):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['sma'])
        return ele_list
    
    def all_ecc_osc(self):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['ecc'])
        return ele_list
    
    def all_inc_osc(self,unit:str='RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['inc'])
        if unit.upper() == 'DEG':
            ele_list = [x*180/np.pi for x in ele_list]
        return ele_list
    
    def all_raan_osc(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['raan'])
        if unit.upper() == 'DEG':
            ele_list = [x*180/np.pi for x in ele_list]
        return ele_list
    
    def all_argp_osc(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            val = dat[5]['argp']
            if val>(2*np.pi-1e-3):
                val = val - 2*np.pi
            ele_list.append(val)
        if unit.upper() == 'DEG':
            ele_list = [x*180/np.pi for x in ele_list]
        return ele_list
    
    def all_M_osc(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            val = dat[5]['M']
            if val>(2*np.pi-1e-3):
                val = val - 2*np.pi
            ele_list.append(val)
        if unit.upper() == 'DEG':
            ele_list = [x*180/np.pi for x in ele_list]
        return ele_list

    def all_n_osc(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['n'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list    
    
    def all_omegatilde_osc(self,unit:str = 'RAD'):
        """ 
        This function returns a list of every dt value contained in the ephemeris object
        """
        ele_list = []
        for dat in self.data:
            ele_list.append(dat[5]['raan'] + dat[5]['argp'])
        if unit.upper() == 'DEG':
            ele_list = [x*np.pi/180 for x in ele_list]
        return ele_list  
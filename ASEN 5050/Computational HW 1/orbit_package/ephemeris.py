import numpy as np
#from orbit_package import orbit
from orbit_package.math_tools import *
import datetime

#orb = orbit()
class Ephemeris():
    def __init__(self,t:list,r_vec: list=None, v_vec:list = None,elements:dict=None,frame:str = None,epoch:datetime.datetime = None,propagation_type: str = None,dt:list=None):
        data = []
        # [t[i],elements[i],r_vec[i],v_vec[i],dt[i],propagation_type]
        if not propagation_type:
            propagation_type = 'NO_PROPAGATION_TYPE'
        for i in range(len(t)):
            sub_list = [t[i]]
            if elements[i]:
                sub_list.append(elements[i])
            else: 
                sub_list.append(None)
            if r_vec[i]:
                sub_list.append(r_vec[i])
            else: 
                sub_list.append(None)
            if v_vec[i]:
                sub_list.append(v_vec[i])
            else: 
                sub_list.append(None)
            if dt[i]:
                sub_list.append(dt[i])
            else: 
                sub_list.append(None)
            sub_list.append(propagation_type)

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
    
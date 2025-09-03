import numpy as np
#from orbit_package import orbit
from math_tools import *
import datetime

#orb = orbit()
class Ephemeris():
    def __init__(self,t:list,elements:dict,frame:str = None,epoch:datetime.datetime = None):
        data = []
        for i in range(len(t)):
            data.append([t[i],elements[i]])
        self.data = data
        self.epoch = epoch
        self.frame = frame
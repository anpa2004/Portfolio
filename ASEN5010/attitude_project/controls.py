import numpy as np

class ControlEvent:
    def __init__(self,t_start,t_end,ref,K,P):
        self.t_start = t_start
        self.t_end = t_end
        self.ref = ref
        self.K = K
        self.P = P

    def is_active(self,t):
        return self.t_start <= t < self.t_end
import numpy as np

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def rotation_matrix(O,w,i):

    R = np.zeros([3,3])

    R[0,0] = cos(w)


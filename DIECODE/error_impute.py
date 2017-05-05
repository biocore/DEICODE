from __future__ import division
import numpy as np
from numpy import linalg as LA
from scipy.linalg import sqrtm, inv
import numpy as np


class error(object):
    
    '''''
    This contains base bench marks for error calculation
    '''''

    def RMSE(org,imputed,mask):
        return (((imputed[mask] - org[mask]) ** 2).mean())**(.5)
    
    def MSE(org,imputed,mask):
        return ((imputed[mask] - org[mask]) ** 2).mean()
    
    def forb(org,imputed,mask):
        return np.linalg.norm(imputed[mask] - org[mask],ord=2)/np.linalg.norm(org[mask],ord=2)
    
    def get_density(data):
        nonzeroscount=np.count_nonzero(data)
        sizel = data.shape
        totalentr=sizel[0]*sizel[1]
        return (nonzeroscount/totalentr*100)

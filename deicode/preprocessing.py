import numpy as np
import pandas as pd
from skbio.stats.composition import closure
from .base import _BaseTransform
np.seterr(all='ignore') #prevent inf errors


class inverse_rclr(_BaseTransform):

    def __init__(self):
        return

    def fit(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)
        self.X_=X_
        self._fit()
        return self

    def _fit(self):
        """ TODO """
        self.X_sp=closure(np.exp(np.array(self.X_.copy()).astype(np.float64)))

    def fit_transform(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)
        self.X_=X_
        self._fit()
        return self.X_sp

class rclr(_BaseTransform):

    def __init__(self,center=True):
        self.center=center
        return

    def fit(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)

        if (X_<0).any():
            raise ValueError('Array Contains Negative Values') 

        self.X_=X_
        self._fit()
        return self

    def _fit(self):
        """ TODO """
        X_log=np.log(np.array(self.X_.copy().astype(float))) #log of all values

        if self.center == True:
        # make masked array
            log_mask = np.array(
                [True] * X_log.shape[0] * X_log.shape[1] 
                ).reshape(X_log.shape)
            log_mask[np.isfinite(X_log)] = False # convert in to zero
            # sum of rows (features)
            m = np.ma.array(X_log, mask=log_mask)
            m=np.subtract(m,m.mean(axis=0))
            m=np.subtract(m.T,m.mean(axis=1).T).T.data
            m[~np.isfinite(X_log)]=np.nan
            self.X_sp=m
        else:
            X_log[~np.isfinite(X_log)]=np.nan
            self.X_sp=X_log
 
    def fit_transform(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)

        if (X_<0).any():
            raise ValueError('Array Contains Negative Values') 
        self.X_=X_
        self._fit()
        return self.X_sp
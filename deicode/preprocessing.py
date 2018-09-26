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
        X_log=np.log(closure(np.array(self.X_.copy().astype(float)))) #log of all values
        log_mask = np.array(
            [True] * X_log.shape[0] * X_log.shape[1] 
            ).reshape(X_log.shape)
        log_mask[np.isfinite(X_log)] = False # convert in to zero
        # sum of rows (features)
        m = np.ma.array(X_log, mask=log_mask)
        gm=m.mean(axis=-1,keepdims=True)
        m = (m - gm).squeeze().data
        m[~np.isfinite(X_log)]=np.nan
        self.X_sp=m
 
    def fit_transform(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)

        if (X_<0).any():
            raise ValueError('Array Contains Negative Values') 
        self.X_=X_
        self._fit()
        return self.X_sp
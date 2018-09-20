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
        X_sp=self.X_
        start=X_sp.copy()
        X_sp=np.log(X_sp) #log of all values
        # make masked array
        logdf_mask = np.array(
            [False] * X_sp.shape[0] * X_sp.shape[1] 
            ).reshape(X_sp.shape)
        logdf_mask[X_sp != -np.inf] = True # convert in to zero
        # sum of rows (features)
        m = np.ma.array(X_sp, mask=logdf_mask)
        beta = m.mean(axis=0)
        X_sp = X_sp - beta
        # sum of columns (samples)
        m = np.ma.array(X_sp, mask=logdf_mask)
        gamma = m.mean(axis=1)
        X_sp = (X_sp.T - gamma.T).T
        X_sp[start==0.0]=np.nan # ensure start values are nan before return
        self.X_sp=X_sp.data.astype(np.float64)
 
    def fit_transform(self,X):
        """ TODO """
        X_=np.array(X.copy()).astype(np.float64)
        self.X_=X_
        self._fit()
        return self.X_sp
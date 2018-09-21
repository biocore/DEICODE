import numpy as np
from numpy import inf
from numpy.matlib import repmat
from numpy.linalg import norm
from scipy.sparse.linalg import svds
from numpy.linalg import matrix_rank
from skbio.stats.composition import (ilr, ilr_inv, clr_inv,
                                     _gram_schmidt_basis)
from deicode._optspace import optspace
from .base import _BaseImpute

class OptSpace(_BaseImpute):

    def __init__(self,rank=None
                ,iteration=5,tol=1e-8
                ,clip_high=True):
        """ TODO """
        self.rank=rank
        self.iteration=iteration
        self.tol=tol
        self.clip_high=clip_high
        return

    def fit(self,X):
        """Fit the model to X_sparse """
        X_sparse=X.copy().astype(np.float64)
        self.X_sparse=X_sparse
        self._fit()
        return self
    
    def _fit(self):
        
        """
            
        Replace all zeros with small non-zero values. Using a collaborative filtering based matrix completion method.
        A replacement for adding only small values to data before performing transforms.
        Also useful for removing sparsity constraints when performing downstream analysis.
        
        ----------
        
        data: array_like a matrix of counts
        rows = Features (i.e. OTUs, metabolites)
        columns = Samples
        
        iteration: float, optional : Default is 40
        The number of convex iterations to optomize the solution
        If iteration is not specified, then the default iteration is 100. Which redcues to a satisfactory error threshold.
        
        
        minval: float, optional : Default is None
        A small number to be used to replace zeros
        If minval is not specified, then the default minval is 1e-3. Worked well in practice with compositional transforms.
        
        clip_high: bool, optional : Default is true
        If maxval is true then values in completed matrix are clipped to be at max the values they exsisted at before imputation
        
        clip_high: float, optional : Default is None
        If maxval is true then values in completed matrix are clipped to be at maxval
        
        Returns
        -------
        numpy.ndarray, np.float64
        A completely dense matrix of counts
        
        
        Raises
        ------
        ValueError
        Raises an error if input is a pandas dataframe and not a numpy array
        `ValueError: Lengths must match to compare`.
        
        
        Notes
        -----
        Assumes a low-rank underlying matrix, this means it performs poorly in gradient like tables. Future high-rank completion methods can overcome this.
        
        
        References
        ----------
        .. [1] Rubinsteyn A, Feldman S. 2016. fancyimpute: Version 0.0.16.
        .. [2] Mazumder R, Hastie T, Tibshirani R. 2010. Spectral Regularization Algorithms for Learning Large Incomplete Matrices. J Mach Learn Res 11:2287–2322.
        .. [3] Pending Publication; Martino and Morton
        
        Examples TODO
        --------
        >>> import numpy as np
        >>> from skbio.stats.impute import complete
        >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
        >>> complete(X)
        array([[ 0.2       ,  0.4       ,  0.4       ,  0.001     ],
        [ 0.23683603,  0.5       ,  0.5       ,  0.00118418]])
        
        """
        
        # make copy for imputation, check type  
        X_sparse=self.X_sparse

        #make rank if none
        if self.rank==None:
            self.rank=matrix_rank(X_sparse)
            if self.rank>=min(X_sparse.shape):
                self.rank=min(X_sparse.shape)-1

        # return solved matrix
        U, s_, V, _  = optspace(X_sparse, r=self.rank, niter=self.iteration, tol=self.tol)
        solution=U.dot(s_).dot(V.T)

        self.solution=solution
        self.feature_weights=V
        self.sample_weights=U
        self.s=s_ 
         
    def fit_transform(self,X):
        """ Returns the solution of fit directly"""
        X_sparse=X.copy().astype(np.float64)
        self.X_sparse=X_sparse
        self._fit()
        return self.solution  
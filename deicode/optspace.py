import numpy as np
from numpy import inf
from numpy.matlib import repmat
from numpy.linalg import norm
from scipy.sparse.linalg import svds
from numpy.linalg import matrix_rank
from scipy.spatial import distance
from skbio.stats.composition import (ilr, ilr_inv, clr_inv,
                                     _gram_schmidt_basis)
from deicode._optspace import optspace
from .base import _BaseImpute

class OptSpace(_BaseImpute):

    def __init__(self,rank=None
                ,iteration=5,tol=1e-5):
        """
            
        OptSpace is a matrix completion algorithm based on a singular value 
        decomposition (SVD) optimized on a local manifold. It has been shown to 
        be quite robust to noise in low rank datasets (1). The objective function that 
        it is trying to optimize over is given by:
            
            min(P|(Y-U*S*V^{T})|_{2}^{2}

        U and V are matrices that are trying to be estimated and S is analogous to a 
        matrix of eigenvalues. Y are the observed values and P is a function such that
        the errors between Y and USV are only computed on the nonzero entries.
        
        Parameters
        ----------
        
        data: numpy.ndarray - a matrix of counts of shape (M,N)
        N = Features (i.e. OTUs, metabolites)
        M = Samples
        
        iteration: float, optional : Default is 5
        The number of convex iterations to optomize the solution
        If iteration is not specified, then the default iteration is 5. Which redcues to a satisfactory error threshold.

        tol: float, optional : Default is 1e-5
        Error reduction break, if the error reduced is less than this value it will return the solution 
        
        Returns
        -------
        U: numpy.ndarray - "Sample Loadings" or the unitary matrix having left singular vectors as columns. Of shape (M,rank)
        s: numpy.ndarray - The singular values, sorted in non-increasing order. Of shape (rank,rank).
        V: numpy.ndarray - "Feature Loadings" or Unitary matrix having right singular vectors as rows. Of shape (N,rank)
        solution: numpy.ndarray - (U*S*V.transpose()) of shape (M,N)
        distance: numpy.ndarray - Distance between each pair of the two collections of inputs. Of shape (M,M)
        
        Raises
        ------
        ValueError

        Raises an error if input is not either pandas dataframe or numpy.ndarray 
            `ValueError: Input data is should be type numpy.ndarray`.

        Raises an error if input shape (M,N) where N>M 
            `ValueError: Data-table contains more samples than features, most likely your data is transposed`.     

        Raises an error if input data does not contain any nans or zeros 
            `ValueError: Data-table contains no missing data in the format np.nan or 0`.

        Raises an error if input data contains infs   
            `ValueError: Data-table contains either np.inf or -np.inf`.

        Raises an error if input data and rank violates min(M,N)<rank   
            `ValueError: The rank must be significantly less than the minimum shape of the input table`.

        References
        ----------
        .. [1] Keshavan RH, Oh S, Montanari A. 2009. Matrix completion from a few entries2009 IEEE International Symposium on Information Theory
        
        Examples 
        --------
        TODO
        
        """

        self.rank=rank
        self.iteration=iteration
        self.tol=tol

        return

    def fit(self,X):
        """
        Fit the model to X_sparse 
        """

        X_sparse=X.copy().astype(np.float64)
        self.X_sparse=X_sparse
        self._fit()
        return self
    
    def _fit(self):
        
        # make copy for imputation, check type  
        X_sparse=self.X_sparse

        if type(X_sparse) is not np.ndarray:
            X_sparse=np.array(X_sparse)
            if type(X_sparse) is not np.ndarray:
                raise ValueError('Input data is should be type numpy.ndarray') 
                
        if X_sparse.shape[0]>X_sparse.shape[1]:
            raise ValueError('Data-table contains more samples than features, most likely your data is transposed') 

        if np.count_nonzero(X_sparse)==0 and np.count_nonzero(~np.isnan(X_sparse))==0: 
            raise ValueError('Data-table contains no missing data in the format np.nan or 0') 

        #make rank if none set a rank 
        if self.rank==None:
            self.rank=matrix_rank(X_sparse)
            if self.rank>=min(X_sparse.shape):
                self.rank=min(X_sparse.shape)-1

        if np.count_nonzero(np.isinf(test_))!=0:
            raise ValueError('Data-table contains either np.inf or -np.inf') 
        
        if self.rank>np.min(X_sparse.shape)
            raise ValueError('The rank must be significantly less than the minimum shape of the input table')

        # return solved matrix
        U, s_, V, _  = optspace(X_sparse, r=self.rank, niter=self.iteration, tol=self.tol)
        solution=U.dot(s_).dot(V.T)

        self.distance=distance.cdist(U,U)
        self.solution=solution
        self.feature_weights=V
        self.sample_weights=U
        self.s=s_ 
         
    def fit_transform(self,X):
        """ 
        Returns the solution of fit directly

        """
        X_sparse=X.copy().astype(np.float64)
        self.X_sparse=X_sparse
        self._fit()
        return self.sample_weights 
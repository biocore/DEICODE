import numpy as np
from numpy import inf
from numpy.matlib import repmat
from numpy.linalg import norm
from scipy.sparse.linalg import svds
from numpy.linalg import matrix_rank
from skbio.stats.composition import (ilr, ilr_inv, clr_inv,
                                     _gram_schmidt_basis)

from DEICODE.opt_space import optspace

class OptSpace(object):

    def __init__(self,data,rank=None,iteration=40,tol=1e-8,minval=None,maxval=None,clip_high=True):
        
        if rank==None:
            rank=matrix_rank(data)
            if rank>=min(data.shape):
                rank=min(data.shape)-1

        self.rank=rank
        self.iteration=iteration
        self.tol=tol
        self.minval=minval
        self.maxval=maxval
        self.clip_high=clip_high
        self.data=data
        
        return

    def complete(self):
        
        
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
        
        Examples
        --------
        >>> import numpy as np
        >>> from skbio.stats.impute import complete
        >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
        >>> complete(X)
        array([[ 0.2       ,  0.4       ,  0.4       ,  0.001     ],
        [ 0.23683603,  0.5       ,  0.5       ,  0.00118418]])
        
        """
        
        otum=self.data.copy().astype(np.float64) # make copy for imputation, check type
        
        # return imputed matrix
        x, s, y, _ = optspace(otum, r=self.rank, niter=self.iteration, tol=self.tol)
        completed=x.dot(s).dot(y.T)
        
        #clip high values
        if self.clip_high==True:
            shape=list(completed.shape)
            for x in range(0, shape[0]):
                for y in range(0, shape[1]):
                    if completed[x, y]>otum[x, y]:
                        completed[x, y]=otum[x, y]

        #clip high values with iter
        if self.maxval is not None and isinstance(self.maxval, (list, tuple, np.ndarray)):
            shape=list(completed.shape)
            for x in range(0, shape[0]):
                for y,maxval_ in zip(range(0, shape[1]),self.maxval):
                    if otum[x, y]==0 and completed[x, y]>maxval_:
                        completed[x, y]=maxval_
    
        #clip high values with float (if both None skip)
        if self.maxval is not None and isinstance(self.maxval, (int, float)):
            shape=list(completed.shape)
            for x in range(0, shape[0]):
                for y in range(0, shape[1]):
                    if otum[x, y]==0 and completed[x, y]>self.maxval:
                        completed[x, y]=self.maxval
        
        #clip values below min
        if self.minval!=None:
            completed[completed<self.minval]=self.minval
        
        return completed

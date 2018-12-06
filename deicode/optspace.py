import numpy as np
from deicode._optspace import optspace
from .base import _BaseImpute
from scipy.spatial import distance
import warnings


class OptSpace(_BaseImpute):

    def __init__(self, rank=2, iteration=5, tol=1e-5):
        """

        OptSpace is a matrix completion algorithm based on a singular value
        decomposition (SVD) optimized on a local manifold. It has been shown to
        be quite robust to noise in low rank datasets (1).
        The objective function that it is trying to optimize
        over is given by:

            min(P|(Y-U*S*V^{T})|_{2}^{2}

        U and V are matrices that are trying to be estimated and S
        is analogous to a matrix of eigenvalues. Y are the
        observed values and P is a function such that
        the errors between Y and USV are only computed
        on the nonzero entries.

        Parameters
        ----------

        X: numpy.ndarray - a rclr preprocessed matrix of shape (M,N)
        N = Features (i.e. OTUs, metabolites)
        M = Samples

        rank: int, optional : Default is 2
        The underlying rank of the default set
        to 2 as the default to prevent overfitting.

        iteration: float, optional : Default is 5
        The number of convex iterations to optomize the solution
        If iteration is not specified, then the default iteration is 5.
        Which redcues to a satisfactory error threshold.

        tol: float, optional : Default is 1e-5
        Error reduction break, if the error reduced is
        less than this value it will return the solution

        Returns
        -------
        U: numpy.ndarray - "Sample Loadings" or the unitary matrix
        having left singular vectors as columns. Of shape (M,rank)

        s: numpy.ndarray - The singular values,
        sorted in non-increasing order. Of shape (rank,rank).

        V: numpy.ndarray - "Feature Loadings" or Unitary matrix
        having right singular vectors as rows. Of shape (N,rank)

        solution: numpy.ndarray - (U*S*V.transpose()) of shape (M,N)

        distance: numpy.ndarray - Distance between each
        pair of the two collections of inputs. Of shape (M,M)

        Raises
        ------
        ValueError

        Raises an error if input is not either dataframe or np.ndarray
            `ValueError: Input data is should be type numpy.ndarray`.

        Raises an error if input data does not contain any nans or zeros
            `ValueError: Data-table contains no missing
            data in the format np.nan or 0`.

        Raises an error if input data contains infs
            `ValueError: Data-table contains either np.inf or -np.inf`.

        Raises an error if input data and rank violates min(M,N)<rank
            `ValueError: The rank must be significantly less than the
            minimum shape of the input table`.

        Raises an error if rank*10> M(Samples)
            `ValueError: There are not sufficient samples to run
            must have rank*10 samples in the table`.

        References
        ----------
        .. [1] Keshavan RH, Oh S, Montanari A. 2009. Matrix completion
                from a few entries (2009_ IEEE International
                Symposium on Information Theory

        Examples
        --------

        >>> from deicode.optspace import OptSpace
        >>> from deicode.preprocessing import rclr
        >>> import numpy as np

        rclr preprocessing

        data is numpy.ndarray
        - a array of counts (samples,features)
        (with shape (M,N) where N>M)

        >>> data=np.array([[3, 3, 0], [0, 4, 2], [3, 0, 1]])
        >>> table_rclr=rclr().fit_transform(data)

        OptSpace (RPCA)

        >>> opt=OptSpace().fit(table_rclr)
        numpy.ndarray - "Sample Loadings"
        >>> U=opt.sample_weights
        numpy.ndarray - "Feature Loadings"
        >>> V=opt.feature_weights
        numpy.ndarray - The singular values
        >>> s=opt.s
        numpy.ndarray - (U*S*V.transpose()) of shape (M,N)
        >>> result=opt.solution

        or

        >>> U,s,V=OptSpace().fit_transform(table_rclr)
        numpy.ndarray - fully dense (no zeros) of shape (M,N)
        >>> result=np.dot(np.dot(U,s),V.T)

        """

        self.rank = rank
        self.iteration = iteration
        self.tol = tol

        return

    def fit(self, X):
        """
        Fit the model to X_sparse
        """

        X_sparse = X.copy().astype(np.float64)
        self.X_sparse = X_sparse
        self._fit()
        return self

    def _fit(self):

        # make copy for imputation, check type
        X_sparse = self.X_sparse

        if not isinstance(X_sparse, np.ndarray):
            X_sparse = np.array(X_sparse)
            if not isinstance(X_sparse, np.ndarray):
                raise ValueError('Input data is should be type numpy.ndarray')

        if (np.count_nonzero(X_sparse) == 0 and
                np.count_nonzero(~np.isnan(X_sparse)) == 0):
            raise ValueError('No missing data in the format np.nan or 0')

        if np.count_nonzero(np.isinf(X_sparse)) != 0:
            raise ValueError('Contains either np.inf or -np.inf')

        if self.rank > np.min(X_sparse.shape):
            raise ValueError('rank must be less than the minimum shape')

        if self.rank * 10 > np.min(X_sparse.shape):
            warnings.warn(
                'Insufficient samples, must have rank*10 samples in the table')

        # return solved matrix
        U, s_, V, _ = optspace(X_sparse, r=self.rank,
                               niter=self.iteration, tol=self.tol)
        solution = U.dot(s_).dot(V.T)

        explained_variance_ = (np.diag(s_) ** 2) / (X_sparse.shape[0] - 1)
        ratio = explained_variance_.sum()
        explained_variance_ratio_ = explained_variance_ / ratio
        self.eigenvalues = np.diag(s_)
        self.explained_variance_ratio = list(explained_variance_ratio_)[::-1]
        self.distance = distance.cdist(U, U)
        self.solution = solution
        self.feature_weights = V
        self.sample_weights = U
        self.s = s_

    def fit_transform(self, X):
        """
        Returns the final SVD of

        U: numpy.ndarray - "Sample Loadings" or the
        unitary matrix having left singular
        vectors as columns. Of shape (M,rank)

        s: numpy.ndarray - The singular values,
        sorted in non-increasing order. Of shape (rank,rank).

        V: numpy.ndarray - "Feature Loadings" or Unitary matrix
        having right singular vectors as rows. Of shape (N,rank)

        """
        X_sparse = X.copy().astype(np.float64)
        self.X_sparse = X_sparse
        self._fit()
        return self.sample_weights, self.s, self.feature_weights

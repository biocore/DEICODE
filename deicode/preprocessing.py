import numpy as np
from skbio.stats.composition import closure
from .base import _BaseTransform
import warnings
np.seterr(all='ignore')
# need to ignore log of zero warning


class rclr(_BaseTransform):

    def __init__(self):
        """

        The rclr procedure first log transform
        the nonzero values before centering the data
        we refer to this preprocessing procedure as
        the robust center log-ratio (rclr) due to its
        ties to the clr transform commonly used in
        compositional data analysis.

        Parameters
        ----------

        X: numpy.ndarray - a table array of all
        positive count data of shape (M,N) containing zeros

        N = Features (i.e. OTUs, metabolites)
        M = Samples

        Returns
        -------

        Raises
        ------
        ValueError

        Raises an error if values in array are negative
            `ValueError: Array Contains Negative Values`.

        Raises an error if Data-table contains either np.inf or -np.inf
            Data-table contains either np.inf or -np.inf

        Raises an error  Data-table contains nans
            Data-table contains nans

        Warning

        RuntimeWarning if there are no zeros in the raw count data
            "Data-table contains no zeros."

        References
        ----------
        -

        Examples
        --------

        >>> from deicode.optspace import OptSpace
        >>> from deicode.preprocessing import rclr
        >>> import numpy as np

        numpy.ndarray - a array of counts (samples,features)
        with shape (M,N) where N>M

        >>> data=np.array([[3, 3, 0], [0, 4, 2], [3, 0, 1]])
        >>> table_rclr=rclr().fit_transform(data)

        """
        return

    def _fit(self):
        """ fits and calc. the rclr  """

        X_ = self.X_.copy().astype(float)

        if (X_ < 0).any():
            raise ValueError('Array Contains Negative Values')

        if np.count_nonzero(np.isinf(X_)) != 0:
            raise ValueError('Data-table contains either np.inf or -np.inf')

        if np.count_nonzero(np.isnan(X_)) != 0:
            raise ValueError('Data-table contains nans')

        if np.count_nonzero(X_) == 0:
            warnings.warn("Data-table contains no zeros.", RuntimeWarning)

        X_log = np.log(closure(np.array(X_)))
        log_mask = np.array(
            [True] * X_log.shape[0] * X_log.shape[1]
        ).reshape(X_log.shape)
        log_mask[np.isfinite(X_log)] = False
        # sum of rows (features)
        m = np.ma.array(X_log, mask=log_mask)
        gm = m.mean(axis=-1, keepdims=True)
        m = (m - gm).squeeze().data
        m[~np.isfinite(X_log)] = np.nan
        self.X_sp = m

    def fit_transform(self, X):
        """ directly returns the rclr transform  """
        X_ = np.array(X.copy()).astype(np.float64)
        self.X_ = X_
        self._fit()
        return self.X_sp

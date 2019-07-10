from abc import abstractmethod


class _BaseImpute(object):

    """Base class for imputation methods.
    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    @abstractmethod
    def fit(self):
        """ Placeholder for fit this
        should be implemented by sub-method"""

    def transform(self):
        """ Apply imputation to X_sparse
        """
        return self.sample_weights

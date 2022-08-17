import unittest
import numpy as np
import numpy.testing as npt
from deicode.preprocessing import rclr
from skbio.stats.composition import clr


class Testpreprocessing(unittest.TestCase):
    def setUp(self):
        self.cdata1 = np.array([[2, 2, 6],
                                [4, 4, 2]])
        self.cdata2 = np.array([[3, 3, 0],
                                [0, 4, 2]])
        self.true2 = np.array([[0.0, 0.0, np.nan],
                               [np.nan, 0.34657359, -0.34657359]])
        self.bad1 = np.array([1, 2, -1])
        self.bad2 = np.array([1, 2, np.inf])
        self.bad3 = np.array([1, 2, np.nan])
        pass

    def test_rclr_dense(self):
        """Test rclr and clr are the same on dense datasets."""
        # test clr works the same if there are no zeros
        cmat = rclr(self.cdata1)
        npt.assert_allclose(cmat, clr(self.cdata1.copy()))

    def test_rclr_sparse(self):
        """Test rclr on sparse data."""
        # test a case with zeros :)
        cmat = rclr(self.cdata2)
        npt.assert_allclose(cmat, self.true2)

    def test_rclr_negative_raises(self):
        """Test rclr ValueError on negative."""
        # test negatives throw value error
        with self.assertRaises(ValueError):
            rclr(self.bad1)

    def test_rclr_inf_raises(self):
        """Test rclr ValueError on undefined."""
        # test undefined throw value error
        with self.assertRaises(ValueError):
            rclr(self.bad2)

    def test_rclr_nan_raises(self):
        """Test rclr ValueError on missing (as nan)."""
        # test nan throw value error
        with self.assertRaises(ValueError):
            rclr(self.bad3)

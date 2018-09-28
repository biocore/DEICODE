from deicode._optspace import (G, F_t, gradF_t, Gp, getoptT, getoptS, optspace)
import numpy as np
from numpy.random import randn, rand
from numpy.linalg import norm
import unittest
import numpy.testing as npt
from gneiss.util import block_diagonal
from scipy.io import loadmat
from skbio.util import get_data_path
from skbio.stats.composition import ilr, ilr_inv
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal


class TestOptspace(unittest.TestCase):
    def setUp(self):
        pass

    def test_G(self):
        X = np.ones((10, 10))
        m0 = 2
        r = 2
        exp = G(X, m0, r)
        self.assertAlmostEqual(exp, 0.644944589179)

    def test_G_z_0(self):
        
        X = np.array([[1,3],[4,1],[2,1]])
        m0 = 2
        r = 2
        exp = G(X, m0, r)
        self.assertAlmostEqual(exp, 2.60980232)

    def test_F_t(self):
        
        X = np.ones((5, 5))
        Y = np.ones((5, 5))
        E = np.zeros((5, 5))
        E[0, 1] = 1
        E[1, 1] = 1
        S = np.eye(5)
        M_E = np.ones((5, 5)) * 6
        M_E[0, 0] = 1
        m0 = 2
        rho = 0.5
        res = F_t(X, Y, S, M_E, E, m0, rho)
        exp = 1
        assert_array_almost_equal(res, exp,decimal=3)

    def test_F_t_random(self):
        
        #random ones and zeros
        np.random.seed(0)
        X = np.ones((5, 5))
        Y = np.ones((5, 5))
        E=np.random.choice([0, 1], size=(5,5))
        S = np.eye(5)
        M_E = np.ones((5, 5)) * 6
        M_E[0, 0] = 1
        m0 = 2
        rho = 0.5
        res = F_t(X, Y, S, M_E, E, m0, rho)
        self.assertAlmostEqual(res, 6.5)

    def test_gradF_t(self):
        X = np.ones((5, 5))
        Y = np.ones((5, 5))
        E = np.zeros((5, 5))
        E[0, 1] = 1
        E[1, 1] = 1
        S = np.eye(5)
        M_E = np.ones((5, 5)) * 6
        M_E[0, 0] = 1
        m0 = 2
        rho = 0.5

        res = gradF_t(X, Y, S, M_E, E, m0, rho)

    def test_Gp(self):
        X = np.ones((5, 5)) * 3
        X[0, 0] = 2
        m0 = 2
        r = 5
        res = Gp(X, m0, r)
        exp = np.array(
            [[1.08731273, 1.6309691, 1.6309691, 1.6309691, 1.6309691],
             [3.57804989, 3.57804989, 3.57804989, 3.57804989, 3.57804989],
             [3.57804989, 3.57804989, 3.57804989, 3.57804989, 3.57804989],
             [3.57804989, 3.57804989, 3.57804989, 3.57804989, 3.57804989],
             [3.57804989, 3.57804989, 3.57804989, 3.57804989, 3.57804989]]
            )

        npt.assert_allclose(exp, res)

    def test_getoptT(self):
        X = np.ones((5, 5))
        Y = np.ones((5, 5))
        E = np.zeros((5, 5))
        E[0, 1] = 1
        E[1, 1] = 1
        S = np.eye(5)
        M_E = np.ones((5, 5)) * 6
        M_E[0, 0] = 1
        m0 = 2
        rho = 0.5
        W, Z = gradF_t(X, Y, S, M_E, E, m0, rho)
        res = getoptT(X, W, Y, Z, S, M_E, E, m0, rho)
        exp = -9.5367431640625e-08
        npt.assert_allclose(exp, res)

    def test_getoptS_small(self):
        # warning : this test must ALWAYS pass
        data = loadmat(get_data_path('small_test.mat'))

        M_E = np.array(data['M_E'].todense())
        E = data['E']

        x = data['x']
        y = data['y']
        res = getoptS(x, y, M_E, E)
        exp = np.array([[ 0.93639499, 0.07644197, -0.02828782],
                        [-0.03960841, 0.60787383, 0.00521257],
                        [ 0.00729038, 0.00785834, 0.67853083]])
        npt.assert_allclose(res, exp, atol=1e-5)

    def test_optspace_original(self):
        M0 = loadmat(get_data_path('large_test.mat'))['M0']
        M_E = loadmat(get_data_path('large_test.mat'))['M_E']

        M0 = M0.astype(np.float)
        M_E = np.array(M_E.todense()).astype(np.float)
        X, S, Y, dist = optspace(M_E, r=3, niter=11, tol=1e-8)
        err = X.dot(S).dot(Y.T) - M0
        n, m = M0.shape

        res = norm(err, 'fro') / np.sqrt(m*n)
        exp = 0.0010701845536
        assert_array_almost_equal(res, exp,decimal=3)

if __name__ == "__main__":
    unittest.main()


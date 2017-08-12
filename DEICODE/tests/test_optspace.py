from DEICODE.optspace import G, F_t, gradF_t, Gp, getoptT, getoptS, optspace
import numpy as np
import unittest
import numpy.testing as npt
from gneiss.util import block_diagonal


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
        # TODO: tests to see when z<1
        pass

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
        self.assertAlmostEqual(res, exp)

    def test_F_t_random(self):
        # TODO: add in test for random zeros
        pass

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
        exp = -0.0125
        npt.assert_allclose(exp, res)

    def test_getoptS(self):
        X = np.ones((5, 5))
        Y = np.ones((5, 5))
        E = np.zeros((5, 5))
        E[0, 1] = 1
        E[1, 1] = 1
        S = np.eye(5)
        M_E = np.ones((5, 5)) * 6
        M_E[0, 0] = 1
        res = getoptS(X, Y, M_E, E)
        exp = np.array([[0.58, 0.58, 0.58, 0.58, 0.58],
                        [0.58, 0.58, 0.58, 0.58, 0.58],
                        [0.58, 0.58, 0.58, 0.58, 0.58],
                        [0.58, 0.58, 0.58, 0.58, 0.58],
                        [0.58, 0.58, 0.58, 0.58, 0.58]])
        npt.assert_allclose(res, exp)

    def test_getoptS_rect(self):
        # TODO: Use rectangular input matrices to test dimensions
        pass

    def test_optspace(self):
        ncols, nrows, nblocks = 10, 10, 2
        np.random.seed(0)
        X_noise = block_diagonal(ncols, nrows, nblocks)
        rng = np.random.RandomState(0)
        spar = 4
        mask = np.random.randint(0, spar, size=X_noise.shape).astype(np.bool)
        rand_zeros = np.random.rand(*X_noise.shape)*0
        X_noise[mask] = rand_zeros[mask]

        res = optspace(X_noise, r=2, niter=10, tol=1e-8)
        print(res)

    def test_optspace_rect(self):
        # TODO: Make the rows and columns different.
        pass

if __name__ == "__main__":
    unittest.main()


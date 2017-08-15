from DEICODE.optspace import G, F_t, gradF_t, Gp, getoptT, getoptS, optspace
import numpy as np
from numpy.random import randn, rand
from numpy.linalg import norm
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

    def test_optspace_no_noise_small(self):
        np.random.seed(0)

        n = 10
        m = 10
        r = 3
        tol = 1e-8

        eps = r*np.log10(n);

        U = randn(n, r)
        V = randn(m, r)
        Sig = np.eye(r)
        M0 = U.dot(Sig).dot(V.T)

        E = 1 - np.ceil( rand(n, m) - eps/np.sqrt(m*n)  )
        M_E = np.multiply(M0, E)
        X, S, Y, dist = optspace(M_E, r=3, niter=10, tol=tol)

        exp_X = np.array([[-8.79865705e-01, 5.80137539e-01, -2.11945417e+00],
                          [2.41254717e+00, -4.15135158e-02, -8.08590887e-02],
                          [9.54706347e-01, 6.53315999e-01, -7.83143452e-01],
                          [-1.26239471e-01, -1.15034917e-01, -6.73104191e-02],
                          [-4.75788895e-04, -1.83763198e-03, 1.77844164e-03],
                          [9.86004154e-01, 1.38993763e+00, 3.02042515e-01],
                          [4.30368032e-03, 1.66220416e-02, -1.60866433e-02],
                          [8.92420611e-01, 3.16881628e-01, 2.81385102e-01],
                          [1.72236780e-01, 1.34451853e+00, -1.60799313e+00],
                          [-8.24468845e-01, 2.31976548e+00, 1.45849936e+00]])

        exp_S = np.array([[0.02717799, -0.00434587, -0.01896239],
                          [-0.23599261, -0.05858438, -0.00468417],
                          [0.20921023, -0.00193656, -0.03990232]])

        exp_Y = np.array([[0.13823077, -0.02226359, 0.16660621],
                          [-0.16332081, -1.83297345, 2.44195478],
                          [0.48775225, 0.12243485, 0.09789257],
                          [0.77352332, -1.52764724, -0.86834177],
                          [1.10757851, 1.45958719, 1.71390276],
                          [-0.42845364, -0.85749063, -0.28152323],
                          [-2.41851741, -0.17623352, 0.19516526],
                          [1.360709, -1.1687075, -0.39642448],
                          [-0.03081843, -0.1683746, 0.18097424],
                          [-0.07726637, -0.00677196, 0.02805902]])

        exp_dist = np.array([0.87510604, 0.87510604, 0.87510604, 0.87510604,
                             0.87510604, 0.87510603, 0.87510603, 0.87510603,
                             0.87510603, 0.87510603, 0.87510603])

        npt.assert_allclose(X, exp_X, atol=1e-5)
        npt.assert_allclose(S, exp_S, atol=1e-5)
        npt.assert_allclose(Y, exp_Y, atol=1e-5)
        npt.assert_allclose(dist, exp_dist, atol=1e-5)

    def test_optspace_noisy_small(self):
        # add some noise
        pass

    def test_optspace_rect(self):
        # TODO: Make the rows and columns different.
        pass

if __name__ == "__main__":
    unittest.main()


import unittest
import numpy as np
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from deicode._optspace import optspace
from skbio.stats.composition import clr
from simulations import build_block_model
from nose.tools import nottest


@nottest
def create_test_table():
    _, test_table = build_block_model(rank=2,
                                      hoced=20,
                                      hsced=20,
                                      spar=2e3,
                                      C_=2e3,
                                      num_samples=50,
                                      num_features=500,
                                      mapping_on=False)

    # the rclr is tested in other places
    # this is just used as input into
    # the OptSpace tests
    test_table = np.array(test_table)
    table_rclr = rclr().fit_transform(test_table)

    return test_table, table_rclr


class TestOptSpace(unittest.TestCase):

    def setUp(self):
        self.test_table, self.test_rclr = create_test_table()
        self.rank = 2
        self.iteration = 5
        self.tol = 1e-5

    def test_OptSpace(self):
        """Tests the basic validity of the
        actual OptSpace() method's output."""

        # run base OptSpace
        opt = OptSpace(rank=self.rank,
                       iteration=self.iteration,
                       tol=self.tol).fit(self.test_rclr)
        U_res, s_res, V_res = OptSpace(
            rank=self.rank,
            iteration=self.iteration,
            tol=self.tol).fit_transform(
            self.test_rclr)
        # use base optspace helper to check
        # that wrapper is not changing outcomes
        U_exp, s_exp, V_exp, dist = optspace(self.test_rclr,
                                             self.rank,
                                             self.iteration,
                                             self.tol)
        # more exact testing of directionally is done
        # in test_method.py. Here we just compare abs
        # see  (c/o @cameronmartino's comment in #29).
        for i in range(self.rank):
            np.testing.assert_array_almost_equal(abs(U_exp[:, i]),
                                                 abs(opt.sample_weights[:, i]))
            np.testing.assert_array_almost_equal(abs(s_exp[:, i]),
                                                 abs(opt.s[:, i]))
            np.testing.assert_array_almost_equal(
                abs(V_exp[:, i]), abs(opt.feature_weights[:, i]))
            np.testing.assert_array_almost_equal(abs(U_exp[:, i]),
                                                 abs(U_res[:, i]))
            np.testing.assert_array_almost_equal(abs(s_exp[:, i]),
                                                 abs(s_res[:, i]))
            np.testing.assert_array_almost_equal(abs(V_exp[:, i]),
                                                 abs(V_res[:, i]))

    def test_OptSpace_rank_raises(self):
        """Tests ValueError for OptSpace() rank."""
        # test rank too low
        try:
            OptSpace(rank=1).fit(self.test_rclr)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")
        # test rank way too high
        try:
            OptSpace(rank=10000).fit(self.test_rclr)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")
        try:
            OptSpace(rank=100).fit(self.test_rclr)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")

    def test_OptSpace_iter_raises(self):
        """Tests ValueError for OptSpace() iteration 0."""
        # test iter too low
        try:
            OptSpace(iteration=0).fit(self.test_rclr)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")

    def test_OptSpace_dense_raises(self):
        """Tests ValueError for OptSpace() on dense and no infs."""
        # test all dense
        try:
            OptSpace(iteration=0).fit(self.test_rclr + 1)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")
        # test wrong clr
        try:
            OptSpace(iteration=0).fit(clr(self.test_table))
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")
        # test wrong datatype
        try:
            OptSpace(iteration=0).fit(1)
        except ValueError:
            pass
        else:
            raise AssertionError("ValueError was not raised")


if __name__ == "__main__":
    unittest.main()

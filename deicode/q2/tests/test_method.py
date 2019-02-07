import unittest
import numpy as np
from biom import Table
from skbio import OrdinationResults
from skbio.stats.distance import DistanceMatrix
from deicode.q2._method import rpca
from simulations import build_block_model


class Testrpca(unittest.TestCase):

    def setUp(self):
        _, test_table = build_block_model(rank=2,
                                          hoced=20,
                                          hsced=20,
                                          spar=2e3,
                                          C_=2e3,
                                          num_samples=50,
                                          num_features=500,
                                          mapping_on=False)

        feat_ids = ['F%d' % i for i in range(test_table.shape[0])]
        samp_ids = ['L%d' % i for i in range(test_table.shape[1])]

        self.test_table = Table(test_table, feat_ids, samp_ids)

    def test_rpca(self):
        ord_test, dist_test = rpca(table=self.test_table)
        self.assertIsInstance(ord_test, OrdinationResults)
        self.assertIsInstance(dist_test, DistanceMatrix)
        self.assertTrue(any(np.isnan(ord_test.features)))
        self.assertTrue(any(np.isnan(ord_test.samples)))


if __name__ == "__main__":
    unittest.main()

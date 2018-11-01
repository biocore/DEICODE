import unittest
import numpy as np
from biom import Table
from skbio import OrdinationResults
from deicode.q2._method import rpca_biplot
from simulations import build_block_model

class Testrpca_biplot(unittest.TestCase):

    def setUp(self):
        _,test_table = build_block_model(rank=2, 
                                 hoced=20, 
                                 hsced=20,
                                 spar=2e3, 
                                 C_=2e3, 
                                 num_samples=200,
                                 num_features=4000,
                                 mapping_on=False)

        feat_ids = ['F%d' % i for i in range(test_table.shape[0])]
        samp_ids = ['L%d' % i for i in range(test_table.shape[1])]

        self.test_table = Table(test_table, feat_ids, samp_ids)

    def test_rpca_biplot(self):
        res = rpca_biplot(table=self.test_table)
        self.assertIsInstance(res, OrdinationResults)
        self.assertTrue(any(np.isnan(res.features)))
        self.assertTrue(any(np.isnan(res.samples)))

if __name__ == "__main__":
    unittest.main()

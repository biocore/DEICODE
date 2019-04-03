from deicode.scripts._standalone_rpca import standalone_rpca
import unittest
import pandas as pd
from click.testing import CliRunner
from skbio import OrdinationResults
from skbio.util import get_data_path
from sklearn.utils.testing import assert_array_almost_equal


class Test_standalone_rpca(unittest.TestCase):
    def setUp(self):
        pass

    def test_standalone_rpca(self):
        in_ = get_data_path('test.biom')
        out_ = '/'.join(in_.split('/')[:-1])
        runner = CliRunner()
        result = runner.invoke(standalone_rpca, ['--in-biom', in_,
                                                 '--output-dir', out_])
        # Read the results
        dist_res = pd.read_csv(get_data_path('distance-matrix.tsv'), sep='\t',
                               index_col=0)
        ord_res = OrdinationResults.read(get_data_path('ordination.txt'))

        # Read the expected results
        dist_exp = pd.read_csv(get_data_path('truth_distance.txt'), sep='\t',
                               index_col=0)
        fea_exp = pd.read_csv(get_data_path('truth_feature.txt'), sep='\t',
                              index_col=0)
        samp_exp = pd.read_csv(get_data_path('truth_sample_trimmed.txt'), sep='\t',
                               index_col=0)

        # Check that the distance matrix matches our expectations
        assert_array_almost_equal(dist_res.values, dist_exp.values)

        # Check that the ordination results match our expectations -- checking
        # each value for both features and samples

        # Row order doesn't matter, but the row names should match
        assert set(ord_res.features.index) == set(fea_exp.index)
        assert set(ord_res.samples.index) == set(samp_exp.index)
        # Column names don't matter, but order does. To enable comparison using
        # pd.testing.assert_frame_equal(), just set the column names to match
        # (once we've verified that there's the same number of columns)
        assert len(ord_res.features.columns) == len(fea_exp.columns)
        assert len(ord_res.samples.columns) == len(samp_exp.columns)
        ord_features_copy = ord_res.features.copy()
        ord_features_copy.columns = fea_exp.columns
        ord_samples_copy = ord_res.samples.copy()
        ord_samples_copy.columns = samp_exp.columns

        # Now, actually check that ordination stuff matches
        pd.testing.assert_frame_equal(ord_features_copy, fea_exp,
                                      check_like=True)
        pd.testing.assert_frame_equal(ord_samples_copy, samp_exp,
                                      check_like=True)

        # Lastly, check that DEICODE's exit code was 0 (indicating success)
        self.assertEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()

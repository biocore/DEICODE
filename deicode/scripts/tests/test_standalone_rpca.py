from deicode.scripts._standalone_rpca import standalone_rpca
import unittest
import pandas as pd
from click.testing import CliRunner
from skbio import OrdinationResults
from skbio.util import get_data_path
from numpy.testing import assert_array_almost_equal


class Test_standalone_rpca(unittest.TestCase):
    def setUp(self):
        pass

    def test_standalone_rpca(self):
        """Checks the distance matrix and biplot produced by DEICODE's script.

           This is more of an "integration test" than a unit test -- the
           details of the algorithm used by the standalone RPCA script are
           checked in more detail in deicode/tests/test_optspace.py, etc.
        """
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
        dist_exp = pd.read_csv(get_data_path('expected-distance-matrix.tsv'),
                               sep='\t', index_col=0)
        ord_exp = OrdinationResults.read(get_data_path(
                                         'expected-ordination.txt'))

        # Check that the distance matrix matches our expectations
        assert_array_almost_equal(dist_res.values, dist_exp.values)

        # Check that the ordination results match our expectations -- checking
        # each value for both features and samples

        # This is just a tuple where the DataFrames within each 2-tuple within
        # it will be compared
        res_exp = ((ord_res.features, ord_exp.features, "feature"),
                   (ord_res.samples, ord_exp.samples, "sample"))

        for (res, exp, aspect) in res_exp:
            # Row order doesn't matter, but the row names should match
            assert set(res.index) == set(exp.index)
            # Column names don't matter, but order does
            assert len(res.columns) == len(exp.columns)
            # Now, we can actually check that the ordination values match.
            # (First, we sort so that the rows in both the result and expected
            # DFs are in the same order, to enable us to use check_names=False
            # when calling pd.testing.assert_series_equal().)
            res_sorted = res.sort_index()
            exp_sorted = exp.sort_index()
            for col_index in range(len(res_sorted.columns)):
                # Extract the n-th column (a pandas Series) from both the
                # result and expected DataFrame, then compare their values.
                res_series = res_sorted.iloc[:, col_index]
                exp_series = exp_sorted.iloc[:, col_index]
                try:
                    # First, just try comparing the two PCs and seeing if their
                    # values are approximately equal.
                    pd.testing.assert_series_equal(res_series, exp_series,
                                                   check_names=False,
                                                   check_less_precise=5)
                except AssertionError:
                    # It's fine for any of the "PC"s (i.e. columns in the
                    # OrdinationResults) to be off by a factor of -1, since
                    # that doesn't really change the interpretation of anything
                    # (c/o @cameronmartino's comment in #29).
                    # To allow for this case to pass the tests, we just try
                    # negating one of the series, and seeing if
                    # that makes them approximately equal.
                    # (If they're *still* not equal, this test will fail.)
                    pd.testing.assert_series_equal(-res_series, exp_series,
                                                   check_names=False,
                                                   check_less_precise=5)
                #print("PC {} for {} ordination matches.".format(col_index,
                #                                                aspect))

        # Lastly, check that DEICODE's exit code was 0 (indicating success)
        self.assertEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()

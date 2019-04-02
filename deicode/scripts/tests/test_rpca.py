from deicode.scripts._rpca import rpca
import unittest
import pandas as pd
from click.testing import CliRunner
from skbio.util import get_data_path
from sklearn.utils.testing import assert_array_almost_equal


class Test_rpca(unittest.TestCase):
    def setUp(self):
        pass

    def test_rpca(self):
        in_ = get_data_path('test.biom')
        out_ = '/'.join(in_.split('/')[:-1])
        runner = CliRunner()
        result = runner.invoke(rpca, ['--in-biom', in_,
                                      '--output-dir', out_])
        dist_res = pd.read_table(get_data_path('distance.txt'), index_col=0)
        fea_res = pd.read_table(get_data_path('feature.txt'), index_col=0)
        samp_res = pd.read_table(get_data_path('sample.txt'), index_col=0)

        dist_exp = pd.read_table(get_data_path('distance.txt'), index_col=0)
        fea_exp = pd.read_table(get_data_path('feature.txt'), index_col=0)
        samp_exp = pd.read_table(get_data_path('sample.txt'), index_col=0)

        assert_array_almost_equal(dist_res.values, dist_exp.values)
        assert_array_almost_equal(fea_res.values, fea_exp.values)
        assert_array_almost_equal(samp_res.values, samp_exp.values)
        self.assertEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()

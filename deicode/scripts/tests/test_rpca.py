from deicode.scripts._rpca import RPCA
import unittest
import pandas as pd
from click.testing import CliRunner
from skbio.util import get_data_path
from sklearn.utils.testing import assert_array_almost_equal


class Test_rpca(unittest.TestCase):
    def setUp(self):
        pass

    def test_RPCA(self):
        in_ = get_data_path('test.biom')
        out_ = '/'.join(in_.split('/')[:-1])
        runner = CliRunner()
        result = runner.invoke(RPCA, ['--in_biom', in_,
                                      '--output_dir', out_])
        dist_res = pd.read_table(get_data_path('distance.txt'), index_col=0)
        fea_res = pd.read_table(get_data_path('feature.txt'), index_col=0)
        samp_res = pd.read_table(get_data_path('sample.txt'), index_col=0)

        dist_exp = pd.read_table(get_data_path('distance.txt'), index_col=0)
        fea_exp = pd.read_table(get_data_path('feature.txt'), index_col=0)
        samp_exp = pd.read_table(get_data_path('sample.txt'), index_col=0)

        assert_array_almost_equal(dist_res.as_matrix(), dist_exp.as_matrix())
        assert_array_almost_equal(fea_res.as_matrix(), fea_exp.as_matrix())
        assert_array_almost_equal(samp_res.as_matrix(), samp_exp.as_matrix())
        self.assertEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()

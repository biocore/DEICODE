from biom import load_table
import pandas as pd
import numpy as np
import os
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
import click


@click.command()
@click.option('--in_biom', help='Input table in biom format.')
@click.option('--output_dir', help='Location of output files.')
@click.option(
    '--min_sample_depth',
    default=500,
    help='Minimum Sample Sequencing Depth Cut Off default=500')
def RPCA(in_biom: str, output_dir: str,
         min_sample_depth: int) -> None:
    """ Runs RPCA with an rclr preprocessing step"""

    table = load_table(in_biom)

    def sample_filter(val, id_, md): return sum(val) > min_sample_depth
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T.drop_duplicates()
    # rclr preprocessing and OptSpace (RPCA)
    opt = OptSpace(rank=3).fit(rclr().fit_transform(table.copy()))
    # Sample Loadings
    sample_loading = pd.DataFrame(opt.sample_weights, index=table.index)
    sample_loading = sample_loading.rename(columns={0: 1, 1: 2, 2: 3})
    # Feature Loadings
    feature_loading = pd.DataFrame(opt.feature_weights, index=table.columns)
    feature_loading = feature_loading.rename(columns={0: 1, 1: 2, 2: 3})
    # make Emperor stype output
    add_ = pd.DataFrame(np.array([[np.nan] * len(sample_loading.columns),
                                  [np.nan] * len(sample_loading.columns),
                                  opt.eigenvalues,
                                  opt.explained_variance_ratio]),
                        index=['',
                               '',
                               'eigvals',
                               '% variation explained'],
                        columns=sample_loading.columns)
    sample_loading = pd.concat([sample_loading, add_], axis=0)
    sample_loading.index.name = 'pc vector number'
    # distance
    pd.DataFrame(opt.distance, table.index,
                 table.index).to_csv(os.path.join(output_dir,
                                                  'distance.txt'), sep='\t')
    # write files to output folder
    feature_loading.to_csv(os.path.join(output_dir, 'feature.txt'), sep='\t')
    sample_loading.to_csv(os.path.join(output_dir, 'sample.txt'), sep='\t')


if __name__ == '__main__':
    RPCA()

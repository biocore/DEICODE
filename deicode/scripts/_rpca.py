from biom import load_table
import pandas as pd
import numpy as np
import skbio
import os
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
import click


@click.command()
@click.option('--in_biom', help='Input table in biom format.')
@click.option('--output_dir', help='Location of output files.')
@click.option('--rank',default=3, help='Rank with witch to run OptSpace. default=3')
@click.option(
    '--min_sample_depth',
    default=500,
    help='Minimum Sample Sequencing Depth Cut Off default=500')


def RPCA(in_biom: str, output_dir: str,
         min_sample_depth: int, rank: int) -> None:
    """ Runs RPCA with an rclr preprocessing step"""

    table = load_table(in_biom)

    def sample_filter(val, id_, md): return sum(val) > min_sample_depth
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T.drop_duplicates()
    # rclr preprocessing and OptSpace (RPCA)
    opt = OptSpace(rank=rank).fit(rclr().fit_transform(table.copy()))
    rename_cols={i-1:'PC'+str(i) for i in range(1,rank+1)}
    
    # Feature Loadings
    feature_loading = pd.DataFrame(opt.feature_weights, index=table.columns)
    feature_loading = feature_loading.rename(columns=rename_cols)

    # Sample Loadings
    sample_loading = pd.DataFrame(opt.sample_weights, index=table.index)
    sample_loading = sample_loading.rename(columns=rename_cols)

    proportion_explained = pd.Series(opt.explained_variance_ratio,
                                     index=list(rename_cols.values()))
    eigvals = pd.Series(opt.eigenvalues,
                        index=list(rename_cols.values()))
    ord_res = skbio.OrdinationResults(
            'PCoA',
            'Principal Coordinate Analysis',
            eigvals,
            sample_loading,
            features=feature_loading,
            proportion_explained=proportion_explained)

    # distance
    pd.DataFrame(opt.distance, table.index,
                 table.index).to_csv(os.path.join(output_dir,
                                                  'Robust_Aitchison_Distance.tsv'), sep='\t')
    # write files to output folder
    ord_res.write(os.path.join(output_dir, 'RPCA_Ordination.txt'))
    return


if __name__ == '__main__':
    RPCA()

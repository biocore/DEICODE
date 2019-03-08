import os
import skbio
import click
import pandas as pd
from biom import load_table
from skbio import OrdinationResults
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr


@click.command()
@click.option('--in_biom', help='Input table in biom format.')
@click.option('--output_dir', help='Location of output files.')
@click.option(
    '--rank',
    default=3,
    help='Rank with witch to run OptSpace. default=3')
@click.option(
    '--min_sample_depth',
    default=500,
    help='Minimum Sample Sequencing Depth Cut Off default=500')
def rpca(in_biom: str, output_dir: str,
         min_sample_depth: int, rank: int) -> None:
    """ Runs RPCA with an rclr preprocessing step"""

    # import table
    table = load_table(in_biom)
    # filter sample to min depth

    def sample_filter(val, id_, md): return sum(val) > min_sample_depth
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T.drop_duplicates()
    # rclr preprocessing and OptSpace (RPCA)
    opt = OptSpace(rank=rank).fit(rclr().fit_transform(table.copy()))
    rename_cols = {i - 1: 'PC' + str(i) for i in range(1, rank + 1)}

    # Feature Loadings
    feature_loading = pd.DataFrame(opt.feature_weights, index=table.columns)
    feature_loading = feature_loading.rename(columns=rename_cols)
    feature_loading.sort_values('PC1', inplace=True, ascending=True)
    feature_loading -= feature_loading.mean(axis=0)

    # Sample Loadings
    sample_loading = pd.DataFrame(opt.sample_weights, index=table.index)
    sample_loading = sample_loading.rename(columns=rename_cols)
    sample_loading -= sample_loading.mean(axis=0)

    proportion_explained = pd.Series(opt.explained_variance_ratio,
                                     index=list(rename_cols.values()))
    eigvals = pd.Series(opt.eigenvalues,
                        index=list(rename_cols.values()))
    # save ordination results
    ord_res = OrdinationResults(
        'PCoA',
        'Principal Coordinate Analysis',
        eigvals.copy(),
        sample_loading.copy(),
        features=feature_loading.copy(),
        proportion_explained=proportion_explained.copy())
    # write files to output folder
    ord_res.write(os.path.join(output_dir, 'RPCA_Ordination.txt'))
    # save distance matrix
    dist_res = skbio.stats.distance.DistanceMatrix(
        opt.distance, ids=sample_loading.index)
    dist_res.write(os.path.join(output_dir, 'RPCA_distance.txt'))
    return


if __name__ == '__main__':
    rpca()

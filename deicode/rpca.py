import biom
import skbio
import pandas as pd
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from deicode._rpca_defaults import (DEFAULT_RANK, DEFAULT_MSC, DEFAULT_MFC,
                                    DEFAULT_ITERATIONS)


def rpca(table: biom.Table,
         rank: int=DEFAULT_RANK,
         min_sample_count: int=DEFAULT_MSC,
         min_feature_count: int=DEFAULT_MFC,
         iterations: int=DEFAULT_ITERATIONS) -> (
         skbio.OrdinationResults,
         skbio.DistanceMatrix):
    """Runs RPCA with an rclr preprocessing step.

       This code will be run by both the standalone and QIIME 2 versions of
       DEICODE.
    """

    # filter sample to min depth
    def sample_filter(val, id_, md): return sum(val) > min_sample_count
    def observation_filter(val, id_, md): return sum(val) > min_feature_count
    table = table.filter(observation_filter, axis='observation')
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T
    if len(table.index) != len(set(table.index)):
        raise ValueError('Data-table contains duplicate indices')
    if len(table.columns) != len(set(table.columns)):
        raise ValueError('Data-table contains duplicate columns')

    # rclr preprocessing and OptSpace (RPCA)
    opt = OptSpace(
        rank=rank,
        iteration=iterations).fit(
        rclr().fit_transform(
            table.copy()))
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

    # % var explained
    proportion_explained = pd.Series(opt.explained_variance_ratio,
                                     index=list(rename_cols.values()))
    # get eigenvalues
    eigvals = pd.Series(opt.eigenvalues,
                        index=list(rename_cols.values()))

    # if the rank is two add PC3 of zeros
    # this is referenced as in issue in
    # <https://github.com/biocore/emperor/commit
    # /a93f029548c421cb0ba365b4294f7a5a6b0209ce>
    # discussed in DEICODE -- PR#29
    if rank == 2:
        feature_loading['PC3'] = [0] * len(feature_loading.index)
        sample_loading['PC3'] = [0] * len(sample_loading.index)
        eigvals.loc['PC3'] = 0
        proportion_explained.loc['PC3'] = 0

    # save ordination results
    short_method_name = 'rpca_biplot'
    long_method_name = '(Robust Aitchison) RPCA Biplot'
    ord_res = skbio.OrdinationResults(
        short_method_name,
        long_method_name,
        eigvals.copy(),
        samples=sample_loading.copy(),
        features=feature_loading.copy(),
        proportion_explained=proportion_explained.copy())
    # save distance matrix
    dist_res = skbio.stats.distance.DistanceMatrix(
        opt.distance, ids=sample_loading.index)

    return ord_res, dist_res

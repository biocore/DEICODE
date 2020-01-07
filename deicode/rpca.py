import biom
import skbio
import numpy as np
import pandas as pd
from deicode.matrix_completion import MatrixCompletion
from deicode.preprocessing import rclr
from deicode._rpca_defaults import (DEFAULT_RANK, DEFAULT_MSC, DEFAULT_MFC,
                                    DEFAULT_ITERATIONS)
from scipy.linalg import svd


def rpca(table: biom.Table,
         n_components: int = DEFAULT_RANK,
         min_sample_count: int = DEFAULT_MSC,
         min_feature_count: int = DEFAULT_MFC,
         max_iterations: int = DEFAULT_ITERATIONS) -> (
        skbio.OrdinationResults,
        skbio.DistanceMatrix):
    """Runs RPCA with an rclr preprocessing step.

       This code will be run by both the standalone and QIIME 2 versions of
       DEICODE.
    """

    # filter sample to min depth
    def sample_filter(val, id_, md): return sum(val) > min_sample_count
    def observation_filter(val, id_, md): return sum(val) > min_feature_count
    # filter and import table
    table = table.filter(observation_filter, axis='observation')
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T
    if len(table.index) != len(set(table.index)):
        raise ValueError('Data-table contains duplicate indices')
    if len(table.columns) != len(set(table.columns)):
        raise ValueError('Data-table contains duplicate columns')

    # rclr preprocessing and OptSpace (RPCA)
    opt = MatrixCompletion(n_components=n_components,
                           max_iterations=max_iterations).fit(rclr(table))
    # get PC column labels for the skbio OrdinationResults
    rename_cols = ['PC' + str(i+1) for i in range(n_components)]
    # get completed matrix for centering
    X = opt.sample_weights @ opt.s @ opt.feature_weights.T
    # center again around zero after completion
    X = X - X.mean(axis=0)
    X = X - X.mean(axis=1).reshape(-1, 1)
    # re-factor the data
    u, s, v = svd(X)
    # only take n-components
    u = u[:, :n_components]
    v = v.T[:, :n_components]
    s = s[:n_components]
    # calc. the new variance using projection
    projection = u * s[:n_components]
    projection_var = np.var(projection, axis=0)
    # get the total approximated variance of the initial table
    tot_var = np.nansum(np.nanvar(rclr(table), axis=1))
    # calculate the fraction of total variance
    p = projection_var / tot_var
    p = p[:n_components]
    # save the loadings
    feature_loading = pd.DataFrame(v, index=table.columns,
                                   columns=rename_cols)
    sample_loading = pd.DataFrame(u, index=table.index,
                                  columns=rename_cols)
    # % var explained
    proportion_explained = pd.Series(p, index=rename_cols)
    # get eigenvalues
    eigvals = pd.Series(s, index=rename_cols)

    # if the n_components is two add PC3 of zeros
    # this is referenced as in issue in
    # <https://github.com/biocore/emperor/commit
    # /a93f029548c421cb0ba365b4294f7a5a6b0209ce>
    # discussed in DEICODE -- PR#29
    if n_components == 2:
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

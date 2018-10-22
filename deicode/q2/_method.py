import biom
import os
import skbio
import numpy as np
import pandas as pd
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr

def rpca_biplot(table: biom.Table, rank: int=3, 
                min_sample_count: int=500, 
                min_feature_count: int=10) -> skbio.OrdinationResults:

    """ Runs RPCA with an rclr preprocessing step"""

    # filter sample to min depth
    def sample_filter(val, id_, md): return sum(val) > min_sample_count
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T.drop_duplicates()
    table = table.T[table.sum()>min_feature_count].T
    table

    # rclr preprocessing and OptSpace (RPCA)
    opt = OptSpace(rank=rank).fit(rclr().fit_transform(table.copy()))
    rename_cols={i-1:'PC'+str(i) for i in range(1,rank+1)}
    
    # Feature Loadings
    feature_loading = pd.DataFrame(opt.feature_weights, index=table.columns)
    feature_loading = feature_loading.rename(columns=rename_cols)
    feature_loading.sort_values('PC1', inplace=True, ascending=True)

    # Sample Loadings
    sample_loading = pd.DataFrame(opt.sample_weights, index=table.index)
    sample_loading = sample_loading.rename(columns=rename_cols)

    proportion_explained = pd.Series(opt.explained_variance_ratio,
                                     index=list(rename_cols.values()))
    eigvals = pd.Series(opt.eigenvalues,
                        index=list(rename_cols.values()))
    # save ordination results 
    short_method_name = 'rpca_biplot'
    long_method_name = '(Robust Aitchison) RPCA Biplot'
    ord_res = skbio.OrdinationResults(short_method_name, long_method_name,
                                      eigvals.copy(),samples=sample_loading.copy(),
                                      features=feature_loading.copy(),
                                      proportion_explained=proportion_explained.copy())

    return ord_res
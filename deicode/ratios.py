import pandas as pd
import numpy as np


def log_ratios(table_tmp, feature_load, sample_load, 
               taxa_tmp=None, axis_sort=0, N_show=3, level='lowest'):
    """"

    Parameters
    ----------

    table_tmp: pandas dataframe - a data table of shape (M,N)
                N = Features (i.e. OTUs, metabolites) - cols
                M = Samples - index
    feature_load: pandas dataframe - a table of shape (P,N)
                M = PC - cols
                N = Features (i.e. OTUs, metabolites) - index
    sample_load: pandas dataframe - a table of shape (P,M)
                M = PC - cols
                N = Samples - index

    taxa_tmp: pandas dataframe - a table of shape (T,M)
                Optional; default use table index
                T = PC - taxa in feature levels
                M = Samples - index
    axis_sort: int - level on loadings to use
                default is 0
    level: str - level of taxon to use
                default is lowest can be any taxa level

    Returns
    -------
    log_ratios: pandas dataframe - log ratios of highly loaded features
    this is in addition to sample loading PC_axis_sort
    metadata catagory and sample IDs

    ratios_: list - The ratios to plot

    Raises
    ------
    ValueError

    Raises an error if axis_sort not in loading axis
        `ValueError: axis_sort must be in loading axis (i.e.) str([list])`.

    Raises an error if input shape (M,N) where N>M
        `ValueError: Data-table contains more samples than features,
        most likely your data is transposed`.

    References
    ----------

    Examples
    --------

    >>> from deicode.ratios import log_ratios
    >>> from deicode.optspace import OptSpace
    >>> from deicode.preprocessing import rclr
    >>> import numpy as np
    >>> import pandas as pd

    rclr preprocessing

    data is numpy.ndarray
    - a array of counts (samples,features)
    (with shape (M,N) where N>M)

    >>> from deicode.ratios import log_ratios
    >>> from deicode.optspace import OptSpace
    >>> from deicode.preprocessing import rclr
    >>> import numpy as np
    >>> import pandas as pd

    >>> data = np.array([[3, 3, 0], [0, 4, 2], [3, 0, 1]])
    >>> meta_data = np.array([['C1', 'C2', 'C3']]).T
    >>> data = pd.DataFrame(data)
    >>> meta_data = pd.DataFrame(meta_data,columns=['Clusters'])
    >>> table_rclr = rclr().fit_transform(data)
    >>> opt = OptSpace().fit(table_rclr)
    >>> sample_weights=pd.DataFrame(opt.sample_weights)
    >>> feature_weights=pd.DataFrame(opt.feature_weights)

    >>> logdf,ratios = log_ratios(data, meta_data,
                                 sample_weights,
                                 feature_weights)

    """

    # make "table"
    if taxa_tmp is None:
        # level has to be lowest but just ID
        level = 'lowest'
        cols_ = ['taxonomy_0', 'taxonomy_1', 'taxonomy_2',
                 'taxonomy_3', 'taxonomy_4', 'taxonomy_5', 'taxonomy_6']
        values_ = [str(c) for c in list(table_tmp.columns)]
        taxa_tmp = pd.DataFrame(np.array([values_] * 7).T,
                                index=table_tmp.columns,
                                columns=cols_)
        taxa_tmp_ = True
    else:
        taxa_tmp=taxa_tmp.astype(str)
        taxa_tmp[taxa_tmp=='nan']=np.nan
        taxa_tmp[taxa_tmp=='None']=np.nan
        taxa_tmp[taxa_tmp=='Unassigned']='__Unclassified'
        taxa_tmp_ = False

    if table_tmp.shape[0] > table_tmp.shape[1]:
        raise ValueError('Data-table contains more samples than features')
    if axis_sort not in feature_load.columns:
        raise ValueError(
            'The axis given to sort is not in the feature rankings')

    # convert unknown taxa to "other"
    taxa_tmp = clean_taxa_table(taxa_tmp, taxa_tmp_)
    # concat features
    feature_taxa = pd.concat([feature_load, taxa_tmp], axis=1).dropna(
        subset=[axis_sort]).sort_values(axis_sort,ascending=False)
    # level groupby
    level_grouping = {level_: feature_taxa.groupby(level_).sum(
    ).sort_values(axis_sort) for level_ in taxa_tmp.columns}
    # get group dicts
    top_otus = bin_level_markers(level_grouping, feature_taxa, level, N_show)
    # get table of ratios
    log_ratios = get_log_ratios(table_tmp, top_otus)
    # add that data back to the dicts
    log_ratios = pd.concat([sample_load, log_ratios], axis=1)
    return log_ratios


def get_taxa(dftmp):
    """ add taxa labels to a certain level given  """
    dftmp.columns = ['kingdom', 'phylum', 'class', 'order',
                     'family', 'genus', 'species']
    return dftmp


def clean_taxa_table(taxa_df, taxa_tmp_):
    """ bin taxa to "lowest" level  and make unclassified as such"""
    mask = np.array([list(taxa_df[l].str.len().values)
                     for l in taxa_df.columns]).T
    mask = pd.DataFrame(mask, taxa_df.index, taxa_df.columns)
    taxa_df[mask < 5] = np.nan
    lowest_assigned = taxa_df.ffill(axis=1).iloc[:, -1]
    if taxa_tmp_:
        taxa_df['lowest'] = [
            'ID:' + str(y) for x,
            y in zip(
                lowest_assigned,
                lowest_assigned.index)]
    else:
        taxa_df['lowest'] = [
            x.split('__')[1] + '(' + x.split('__')[0] + ')_{ID:' + str(y) + '}'
            for x, y in zip(
                lowest_assigned,
                lowest_assigned.index)]
    taxa_df[mask < 5] = '__Unclassified'

    return taxa_df


def bin_level_markers(level_gp, taxmatch, level, N_show):
    """otus by taxa level given for N x_n,y_n by feature ranking"""
    ratios = {}
    for i in range(N_show):
        top_ = level_gp[level].iloc[[i]].index[0]
        bottom_ = level_gp[level].iloc[[-(i + 1)]].index[0]
        top_l = list(taxmatch[taxmatch[level].isin([top_])].index)
        bottom_l = list(taxmatch[taxmatch[level].isin([bottom_])].index)       
        ratios[(top_,bottom_)]= (top_l,bottom_l)

    return ratios


def get_log_ratios(tabledf, topd):
    """ get log ratios for observed taxa given by bin_level_markers"""
    log_ratios = []
    for (x_i,y_j),(x_i_features,y_j_features) in topd.items():
        if not isinstance(x_i_features, (list,)):
            x_i_features=[x_i_features]
        if not isinstance(y_j_features, (list,)):
            y_j_features=[y_j_features]
        p_x_i_y_j = tabledf.loc[:, list(x_i_features)+list(y_j_features)]
        x_i_vector = tabledf.loc[:, list(x_i_features)][p_x_i_y_j.T.min() > 0]
        y_j_vector = tabledf.loc[:, list(y_j_features)][p_x_i_y_j.T.min() > 0]
        tmp_ratio = (np.log(x_i_vector).mean(axis=1) -
                     np.log(y_j_vector).mean(axis=1))
        col_ = [
            'log(\dfrac{' +
            x_i.replace(
                '__',
                '') +
            '}{' +
            y_j.replace(
                '__',
                '') +
            '})']
        log_ratios.append(pd.DataFrame(tmp_ratio, columns=col_))
    return pd.concat(log_ratios, axis=1)
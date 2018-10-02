from __future__ import division
# utils
import pandas as pd
import numpy as np
from collections import Counter
# blocks
from scipy.stats import norm
from numpy.random import poisson, lognormal
from skbio.stats.composition import closure
from scipy.special import kl_div
# minimize model perams
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
# Set random state
rand = np.random.RandomState(42)


def get_taxa(dftmp):
    """ TODO """
    dftmp.columns = ['kingdom', 'phylum', 'class', 'order',
                     'family', 'genus', 'species']
    dftmp['taxonomy'] = dftmp[dftmp.columns].apply(
        lambda x: ';'.join(x), axis=1)
    return dftmp


def get_enriched_labels(feature_plot, cutoff=0):
    """ TODO """
    # make copies
    feature_plot = feature_plot.copy()

    # round to avoid too many otus
    feature_plot = np.around(feature_plot, 1)
    # get N otus to view
    N_otu = np.min(
        (dict(
            Counter(
                feature_plot[0] < -
                cutoff))[True],
         dict(
            Counter(
                feature_plot[0] > cutoff))[True]))
    feature_plot_low = feature_plot[:N_otu].index
    feature_plot_high = feature_plot[feature_plot.shape[0] - N_otu:].index

    return list(feature_plot_low) + list(feature_plot_high)


def get_mean_abundance(plot_table, plot_map, features_enrich, plot_col):
    """ TODO """
    plot_table = plot_table.copy()
    plot_map = plot_map.copy()
    plot_table = pd.DataFrame(
        closure(plot_table) * 100,
        plot_table.index,
        plot_table.columns).T.copy()
    plot_table.columns = plot_map[plot_col]
    plot_table = plot_table.T.groupby(plot_table.columns).mean()[
        features_enrich].T
    return plot_table


def get_lowest_level(taxonomy, features_en, default_highest=-3):
    """ TODO """
    lowest_level = {}
    taxonomy_enriched = taxonomy.loc[features_en, :].copy()
    for level_ in taxonomy_enriched.columns[:default_highest]:
        tax_tmp = taxonomy_enriched[taxonomy_enriched[level_].apply(
            lambda x: len(x) > 3)]
        for ind_ in taxonomy_enriched[taxonomy_enriched[level_].apply(
                lambda x: len(x) > 3)].index:
            name_ = tax_tmp.loc[ind_, level_].split('__')
            name_ = name_[1] + ' (' + name_[0] + ')'
            lowest_level[ind_] = name_
    return lowest_level


def mean_KL(a, b):
    """
    Returns the KL divergence between matrix A and B

    Parameters
    ----------

    a : dataframe

    b : dataframe


    Returns
    -------
    KL Divergence as list over samples
    """

    kl = []
    a = np.array(a.copy())
    b = np.array(b.copy())
    for i in range(a.shape[0]):
        kl += [kl_div(a[i], b[i])]
    return np.mean(kl)


def Homoscedastic(X_noise, intensity):
    """ uniform normally dist. noise """
    X_noise = np.array(X_noise)
    err = intensity * np.ones_like(X_noise.copy())
    X_noise = rand.normal(X_noise.copy(), err)

    return X_noise


def Heteroscedastic(X_noise, intensity):
    """ non-uniform normally dist. noise """
    err = intensity * np.ones_like(X_noise)
    i = rand.randint(0, err.shape[0], 5000)
    j = rand.randint(0, err.shape[1], 5000)
    err[i, j] = intensity
    X_noise = abs(rand.normal(X_noise, err))

    return X_noise


def Subsample(X_noise, spar, num_samples):
    """ yij ~ PLN( lambda_{ij}, /phi ) """
    # subsample
    mu = spar * closure(X_noise.T).T
    X_noise = np.vstack([poisson(lognormal(np.log(mu[:, i]), 1))
                         for i in range(num_samples)]).T
    # add sparsity

    return X_noise


def block_diagonal_gaus(
        ncols,
        nrows,
        nblocks,
        overlap=0,
        minval=0,
        maxval=1.0):
    """
    Generate block diagonal with Gaussian distributed values within blocks.

    Parameters
    ----------

    ncol : int
        Number of columns

    nrows : int
        Number of rows

    nblocks : int
        Number of blocks, mucst be greater than one

    overlap : int
        The Number of overlapping columns (Default = 0)

    minval : int
        The min value output of the table (Default = 0)

    maxval : int
        The max value output of the table (Default = 1)


    Returns
    -------
    np.array
        Table with a block diagonal where the rows represent samples
        and the columns represent features.  The values within the blocks
        are gaussian distributed between 0 and 1.
    Note
    ----
    The number of blocks specified by `nblocks` needs to be greater than 1.

    """

    if nblocks <= 1:
        raise ValueError('`nblocks` needs to be greater than 1.')
    mat = np.zeros((nrows, ncols))
    gradient = np.linspace(0, 10, nrows)
    mu = np.linspace(0, 10, ncols)
    sigma = 1
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    mat = np.vstack(xs).T

    block_cols = ncols // nblocks
    block_rows = nrows // nblocks
    for b in range(nblocks - 1):

        gradient = np.linspace(5, 5, block_rows)  # samples (bock_rows)
        # features (block_cols+overlap)
        mu = np.linspace(0, 10, block_cols + overlap)
        sigma = 2.0
        xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
              for i in range(len(mu))]

        B = np.vstack(xs).T * maxval
        lower_row = block_rows * b
        upper_row = min(block_rows * (b + 1), nrows)
        lower_col = block_cols * b
        upper_col = min(block_cols * (b + 1), ncols)

        if b == 0:
            mat[lower_row:upper_row, lower_col:int(upper_col + overlap)] = B
        else:
            ov_tmp = int(overlap / 2)
            if (B.shape) == (mat[lower_row:upper_row,
                             int(lower_col-ov_tmp):int(upper_col+ov_tmp+1)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col + ov_tmp + 1)] = B
            elif (B.shape) == (mat[lower_row:upper_row,
                               int(lower_col - ov_tmp):int(upper_col + ov_tmp)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col + ov_tmp)] = B
            elif (B.shape) == (mat[lower_row:upper_row, 
                               int(lower_col - ov_tmp):int(upper_col + ov_tmp - 1)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col +ov_tmp - 1)] = B

    upper_col = int(upper_col - overlap)
    # Make last block fill in the remainder
    gradient = np.linspace(5, 5, nrows - upper_row)
    mu = np.linspace(0, 10, ncols - upper_col)
    sigma = 4
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    B = np.vstack(xs).T * maxval

    mat[upper_row:, upper_col:] = B

    return mat


def build_block_model(
        rank,
        hoced,
        hsced,
        spar,
        C_,
        num_samples,
        num_features,
        overlap=0,
        mapping_on=True):
    """
    Generates hetero and homo scedastic noise on base truth block diagonal with Gaussian distributed values within blocks.

    Parameters
    ----------

    rank : int
        Number of blocks


    hoced : int
        Amount of homoscedastic noise

    hsced : int
        Amount of heteroscedastic noise

    inten : int
        Intensity of the noise

    spar : int
        Level of sparsity

    C_ : int
        Intensity of real values

    num_features : int
        Number of rows

    num_samples : int
        Number of columns

    overlap : int
        The Number of overlapping columns (Default = 0)

    mapping_on : bool
        if true will return pandas dataframe mock mapping file by block


    Returns
    -------
    Pandas Dataframes
    Table with a block diagonal where the rows represent samples
    and the columns represent features.  The values within the blocks
    are gaussian.

    Note
    ----
    The number of blocks specified by `nblocks` needs to be greater than 1.

    """

    # make a mock OTU table
    X_true = block_diagonal_gaus(
        num_samples,
        num_features,
        rank,
        overlap,
        minval=.01,
        maxval=C_)
    if mapping_on:
        # make a mock mapping data
        mappning_ = pd.DataFrame(np.array([['Cluster %s' %
                                            str(x)] *
                                           int(num_samples /
                                               rank) for x in range(1, rank +
                                                                    1)]).flatten(), columns=['example'], index=['sample_' +
                                                                                                                str(x) for x in range(0, num_samples -
                                                                                                                                      2)])

    X_noise = X_true.copy()
    X_noise = np.array(X_noise)
    # add Homoscedastic noise
    X_noise = Homoscedastic(X_noise, hoced)
    # add Heteroscedastic noise
    X_noise = Heteroscedastic(X_noise, hsced)
    # Induce low-density into the matrix
    X_noise = Subsample(X_noise, spar, num_samples)

    # return the base truth and noisy data
    if mapping_on:
        return X_true, X_noise, mappning_
    else:
        return X_true, X_noise


def minimize_model(x0, bnds, X_true_):
    """ TODO """
    def add_noise_min(x_o, X_true=X_true_):
        """ TODO """
        rank = 3
        hoced = x_o[1]
        hsced = x_o[2]
        spar = x_o[3]
        C_ = x_o[4]
        overlap_ = x_o[5]

        num_samples = int(X_true.shape[0])
        num_features = int(X_true.shape[1])

        # make a mock OTU table
        X_true_ = block_diagonal_gaus(
            num_samples,
            num_features,
            rank,
            overlap=overlap_,
            minval=.01,
            maxval=C_)

        X_noise = X_true.copy()
        X_noise = np.array(X_noise)
        # add Homoscedastic noise
        X_noise = Homoscedastic(X_noise, hoced)
        # add Heteroscedastic noise
        X_noise = Heteroscedastic(X_noise, hsced)
        # Induce low-density into the matrix
        X_noise = Subsample(X_noise, spar, num_samples)

        # return the noisy data
        return mean_squared_error(X_true.T, X_noise)

    model_fit = minimize(add_noise_min, x0, bounds=bnds)
    return model_fit


def build_grad_model(hoced, hsced, spar, sigma, C_, num_samples, num_features):
    """ TODO """

    gradient = np.linspace(0, 10, num_samples)
    mu = np.linspace(0, 10, num_features)
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    table = np.vstack(xs).T
    table = pd.DataFrame(table)
    table.index = table.index.astype(np.str)
    table.columns = table.columns.astype(np.str)
    table = table * C_
    X_true = table.values.T.copy()

    X_noise = X_true.copy()
    X_noise = np.array(X_noise)
    # add Homoscedastic noise
    X_noise = Homoscedastic(X_noise, hoced)
    # add Heteroscedastic noise
    X_noise = Heteroscedastic(X_noise, hsced)
    # Induce low-density into the matrix
    X_noise = Subsample(X_noise, spar, num_samples)

    # make a mock mapping data
    mappning_ = pd.DataFrame(np.array([x for x in range(0, num_samples)]), columns=[
                             'Gradient'], index=['sample_' + str(x) for x in range(0, num_samples)])
    mappning_ = mappning_.apply(pd.to_numeric, errors='ignore')

    # return the base truth and noisy data
    return X_true, X_noise, mappning_


def minimize_model_grad(x0, bnds, X_true_):
    """ TODO """

    def add_noise_min_grad(x_o, X_true=X_true_):
        """ TODO """
        # build model and minmize kl-div

        hoced = x_o[0]
        hsced = x_o[1]
        spar = x_o[2]
        sigma = x_o[3]
        C_ = x_o[4]

        num_samples = int(X_true.shape[0])
        num_features = int(X_true.shape[1])

        gradient = np.linspace(0, 10, num_samples)
        mu = np.linspace(0, 10, num_features)
        xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
              for i in range(len(mu))]
        table = np.vstack(xs).T
        table = pd.DataFrame(table)
        table.index = table.index.astype(np.str)
        table.columns = table.columns.astype(np.str)
        table = table * C_

        X_noise = table.values.T.copy()
        X_noise = np.array(X_noise)
        # add Homoscedastic noise
        X_noise = Homoscedastic(X_noise, hoced)
        # add Heteroscedastic noise
        X_noise = Heteroscedastic(X_noise, hsced)
        # Induce low-density into the matrix
        X_noise = Subsample(X_noise, spar, num_samples)

        # return the noisy data
        return mean_squared_error(X_true.T, X_noise)

    model_fit = minimize(add_noise_min_grad, x0, bounds=bnds)

    return model_fit

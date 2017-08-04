from skbio.stats.composition import closure
from gneiss.composition import ilr_transform
from gneiss.regression import ols
from gneiss.cluster import correlation_linkage
from fancyimpute import SoftImpute
from sklearn.model_selection import KFold
from gneiss.util import _type_cast_to_float
import pandas as pd
import skbio
import numpy as np
from patsy import dmatrix


def ols_regression(table: pd.DataFrame, tree: skbio.TreeNode,
                   metadata: pd.DataFrame, formula: str,
                   complete_matrix: bool) -> None:
    """
    Parameters
    ----------
    table : pd.DataFrame
       Input table where samples are rows, features are columns.
    tree : skbio.TreeNode
       Tree object used to perform the ilr transform.
    formula : str
       Patsy formula to specify regression model.

    Returns
    -------
    res: OLSModel
       Optimal fitted regression object for performing the OLS.
    params: pd.Series
       Optimal parameters from hyper parameter estimation.
    cv: pd.DataFrame
       Cross validation results.
    """
    metadata = _type_cast_to_float(metadata)
    num_folds = 5
    if complete_matrix == True:
        res, beta, alpha, mse_score = optimize_complete(table, tree,
                                                        formula, metadata)
        params = pd.Series({'beta': beta, 'alpha': alpha, 'mse_score': mse_score})
    else:
        balances = ilr_transform(table, tree=tree)
        res = ols(table=balances, metadata=metadata,
                  formula=formula)
        res.fit()
        params = pd.Series()

    # perform cross validation
    if complete_matrix == False:
        cv = res.kfold(num_folds)
    else:
        cv = pd.DataFrame(index=['fold_%d' % i for i in range(num_folds)],
                          columns=['model_mse', 'Rsquared', 'pred_mse'],
                          dtype=np.float64)
        kf = KFold(n_splits=num_folds, shuffle=True)
        for i, (train, test) in enumerate(kf.split(np.arange(len(table.index)))):
            Y_train, Y_test = table.iloc[train], table.iloc[test]
            X_train, X_test = metadata.iloc[train], metadata.iloc[test]
            res_, beta, alpha, mse_score = optimize_complete(
                Y_train, tree, formula, X_train)

            X_train = dmatrix(formula, X_train, return_type='dataframe')
            X_test = dmatrix(formula, X_test, return_type='dataframe')
            # model error
            p = res_.predict(X=X_train).values
            y = pd.DataFrame(complete(closure(Y_train.as_matrix().copy()),
                                      int(beta), alpha),
                             columns=Y_train.columns, index=Y_train.index)
            r = ilr_transform(y, tree=tree)

            model_resid = ((p - r)**2)
            model_mse = np.mean(model_resid.sum(axis=0))

            cv.loc['fold_%d' % i, 'model_mse'] = model_mse
            cv.loc['fold_%d' % i, 'Rsquared'] = res_.r2

            # prediction error
            p = res_.predict(X=X_test).values
            y = pd.DataFrame(complete(closure(Y_test.as_matrix().copy()),
                                      int(beta), alpha),
                             columns=Y_test.columns, index=Y_test.index)
            r = ilr_transform(y, tree=tree)

            pred_resid = ((p - r)**2)
            pred_mse = np.mean(pred_resid.sum(axis=0))
            cv.loc['fold_%d' % i, 'pred_mse'] = pred_mse

    return res, params, cv


def optimize_complete(df_, tree_, formula_, metadata_, min_alpha=1e-3, max_alpha=1,
                      alpha_iter=2000, min_beta=1, max_beta=1000,
                      beta_iter=10, iter_max=200, in_thresh=.001):
    """
    Takes OTU table, formula and metadata and creates a
    hill-climbing optimization based on minimizing ols mse.
    First optimizes alpha at a low beta and then proceeds to
    optimize computationaly costly beta operations.

    ----------

    df_: Dataframe matrix of counts
    rows = Features (i.e. OTUs, metabolites)
    columns = Samples

    formula_: str, ols formula linked to metadata

    metadata_: Dataframe
         metadata matched to df_
    rows = Samples
    columns = Mapping data titles

    min_alpha: numpy.float
        minimum search value for alpha (min completion values), default: 1e-3
    max_alpha: numpy.float
        max search value for alpha (min completion values)
        default: 1 (should not exceed 1)
    alpha_iter: numpy.int
        number of values to search between min and max alpha, default: 50
    min_beta: numpy.int
        minimum search value for beta (min number of completion iterations)
        default: 5
    max_beta: numpy.int
        max search value for beta (max number of completion iterations)
        default: 1000
    beta_iter: numpy.int
        number of values to search between min and max beta, default: 10
    iter_max: numpy.int
        maximum number of climbing iterations to search for beta and alpha
        default: 200
    in_thresh: numpy.float
        threshold for difference between previous old mse and current mse
        before terminating climb, default: .001

    Returns
    -------
    Dataframe, completed table of counts
    List, of alpha input searched
    List, of beta input searched
    np.float, optimized mse score


    Raises
    ------
    ValueError
    Raises an error if max alpha > 1:
    `ValueError: max alpha must be less than or equal to 1`.

    ValueError
    Raises an error if alpha_iter or beta_iter exceed iter_max:
    `ValueError: alpha_iter or beta_iter exceed iter_max`.
    """
    alpha_search=np.linspace(min_alpha, max_alpha, alpha_iter)
    beta_search=np.linspace(min_beta, max_beta, beta_iter, dtype=int)
    return climb_(df_, formula_, metadata_, tree_,
                  alpha_search, beta_search, iter_max, in_thresh)


def complete(data, iteration, minval):

    """

    Replace all zeros with small non-zero values. Using a collaborative filtering
    based matrix completion method.  A replacement for adding only small values
    to data before performing transforms. Also useful for removing sparsity
    constraints when performing downstream analysis.

    ----------

    data: array_like a matrix of counts
    rows = Features (i.e. OTUs,  metabolites)
    columns = Samples

    iteration: float,  optional
    The number of convex iterations to optomize the solution
    If iteration is not specified,  then the default iteration is 100.
    Which redcues to a satisfactory error threshold.


    minval: float,  optional
    A small number to be used to replace zeros
    If minval is not specified,  then the default minval is 1e-3.
    Worked well in practice with compositional transforms.



    Returns
    -------
    numpy.ndarray,  np.float64
    A completely dense matrix of counts


    Raises
    ------
    ValueError
    Raises an error if input is a pandas dataframe and not a numpy array
    `ValueError: Lengths must match to compare`.


    Notes
    -----
    Assumes a low-rank underlying matrix,  this means it performs poorly in
    gradient like tables. Future high-rank completion methods can overcome this.


    References
    ----------
    .. [1] Rubinsteyn A,  Feldman S. 2016. fancyimpute: Version 0.0.16.
    .. [2] Mazumder R,  Hastie T,  Tibshirani R. 2010. Spectral Regularization
           Algorithms for Learning Large Incomplete Matrices. J Mach Learn
           Res 11:2287–2322.
    .. [3] Pending Publication; Martino and Morton

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.impute import complete
    >>> X = np.array([[.2, .4, .4,  0], [0, .5, .5, 0]])
    >>> complete(X)
    array([[ 0.2       ,   0.4       ,   0.4       ,   0.001     ],
    [ 0.23683603,   0.5       ,   0.5       ,   0.00118418]])

    """

    otum=data.copy().astype(np.float64) # make copy for imputation,  check type
    otum[otum == 0] = np.nan # make all previously zero values unknown (i.e. Nan)
    #return imputed matrix through Fancy Impute's Soft Impute function
    return SoftImpute(max_rank=min(otum.shape), max_iters=iteration,
                      convergence_threshold=0.00001,
                      min_value=minval, max_value=(np.amax(otum)),
                      verbose=False).complete(otum)


def climb_(df_, formula_, metadata_, tree_,
           alpha, beta, iter_max=200, in_thresh=.001):

    #Find ALPHA
    iter_=1
    res_complete, score_=ols_score(df_, formula_, metadata_, alpha[0],  10,
                                   tree_=tree_)
    for alpha_ in alpha[1:]:
        res_complete, score_tmp=ols_score(df_, formula_, metadata_, alpha_,  10,
                                          tree_=tree_)
        thresh_=score_-score_tmp
        score_=score_tmp
        iter_+=1
        if iter_ > iter_max or (thresh_ <= in_thresh):
            break

    #Use alpha to find BETA
    iter_=1
    res_complete, score_=ols_score(df_, formula_, metadata_, alpha_,  beta[0],
                                   tree_=tree_)
    for beta_ in beta[1:]:
        res_complete, score_tmp=ols_score(df_, formula_, metadata_, alpha_,  beta_,
                                          tree_=tree_)
        thresh_ = score_ - score_tmp
        score_ = score_tmp
        iter_+=1
        if iter_ > iter_max or (thresh_ <= in_thresh):
            break

    return res_complete, beta_, alpha_, score_


def ols_score(df_, formula_, metadata_, alpha, beta,
              build_tree=False, tree_=None):
    """

    Takes OTU table, formula, metadata, alpha and beta values and returns
    a mse score.

    Parameters
    ----------
    df_: Dataframe matrix of counts
    rows = Features (i.e. OTUs, metabolites)
    columns = Samples

    formula_: str, ols formula linked to metadata

    metadata_: Dataframe, metadata matched to df_
    rows = Samples
    columns = Mapping data titles

    alpha: numpy.float, minimum value for completed matrix

    beta: numpy.int, number of iterations for completion optimization

    Returns
    -------
    numpy.float, ols regression mse sum.

    """
    completed=pd.DataFrame(complete(closure(df_.as_matrix().copy()),
                                    int(beta), alpha),
                           columns=df_.columns, index=df_.index)

    balances = ilr_transform(completed, tree_)
    res_complete = ols(formula_, balances, metadata_)
    res_complete.fit()
    return res_complete, res_complete.mse.sum()

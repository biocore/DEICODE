#utils
from __future__ import division
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
#PCoA
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
#blocks
from scipy.stats import norm
from numpy.random import poisson, lognormal
from skbio.stats.composition import closure
from scipy import stats
from scipy.special import kl_div
from sklearn.metrics import mutual_info_score

#minimize model perams
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
# Set random state
rand = np.random.RandomState(42)

def rclr(logdf:'Input dataframe or numpy array of all features in rows and all samples in columns')->'Output is robust clr':
    np.seterr(all='ignore')
    if isinstance(logdf, pd.DataFrame):
        logdf=logdf.copy().as_matrix() # covert df 
        start=logdf.copy()
        logdf=np.log(logdf) #log of all values
    else:
        start=logdf.copy()
        logdf=np.log(logdf) #log of all values

    logdf_mask = np.array(
        [False] * logdf.shape[0] * logdf.shape[1] 
        ).reshape(logdf.shape)
    logdf_mask[logdf_mask != -np.inf] = True # convert in to zero
    # sum of rows (features)
    m = np.ma.array(logdf, mask=logdf_mask)
    beta = m.mean(axis=0)
    logdf = logdf - beta

    # sum of columns (samples)
    m = np.ma.array(logdf, mask=logdf_mask)
    gamma = m.mean(axis=1)
    rclr = (logdf.T - gamma.T).T

    rclr[start==0.0]=np.nan # ensure start values are nan before return
    
    return rclr.data.astype(np.float64)

def Homoscedastic(X_noise,intensity):
    
    X_noise=np.array(X_noise)
    err = intensity * np.ones_like(X_noise.copy())
    X_noise = rand.normal(X_noise.copy(), err)
    
    return X_noise


def Heteroscedastic(X_noise,intensity):
    
    err = intensity * np.ones_like(X_noise)
    i = rand.randint(0, err.shape[0], 5000)
    j = rand.randint(0, err.shape[1], 5000)
    err[i, j] = intensity
    X_noise = abs(rand.normal(X_noise, err))
    
    return X_noise

def Subsample(X_noise,spar,num_samples):
    
    # subsample
    mu=spar*closure(X_noise.T).T
    X_noise=np.vstack([poisson(lognormal(np.log(mu[:,i]),1)) for i in range(num_samples)]).T
    # add sparsity
    
    
    return X_noise


def block_diagonal_gaus(ncols,nrows,nblocks,overlap=0,minval=0,maxval=1.0):
    
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
    
    
    rank = 2**3 - 1
    gradient = np.linspace(0, 10, nrows)
    mu = np.linspace(0, 10, ncols)
    sigma = 1
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    mat = np.vstack(xs).T
    
    
    block_cols = ncols // nblocks
    block_rows = nrows // nblocks
    for b in range(nblocks-1):
    
        gradient = np.linspace(5, 5, block_rows) # samples (bock_rows)
        mu = np.linspace(0, 10, block_cols+overlap) # features (block_cols+overlap)
        sigma = 2.0
        xs = [norm.pdf(gradient, loc=mu[i], scale=sigma) for i in range(len(mu))]
        
        B = np.vstack(xs).T*maxval
                
        #B = np.random.uniform(minval,maxval,size=(block_rows, block_cols))
        lower_row = block_rows * b
        upper_row = min(block_rows*(b+1), nrows)
        lower_col = block_cols * b
        upper_col = min(block_cols*(b+1), ncols)
        

        
        if b==0:
            mat[lower_row:upper_row, lower_col:int(upper_col+overlap)] = B
        else:
            if (B.shape)==(mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2)+1)].shape):    
                mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2)+1)] = B
            elif (B.shape)==(mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2))].shape):    
                mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2))] = B
            elif (B.shape)==(mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2)-1)].shape):    
                mat[lower_row:upper_row, int(lower_col-int(overlap/2)):int(upper_col+int(overlap/2)-1)] = B

    upper_col=int(upper_col-overlap)
    # Make last block fill in the remainder
    gradient = np.linspace(5, 5, nrows-upper_row)
    mu = np.linspace(0, 10, ncols-upper_col)
    sigma = 4
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    B = np.vstack(xs).T*maxval
    
    
    mat[upper_row:, upper_col:] = B
    
    
    return mat

def plot_pcoa(otupcoadfs,mapping,catvis,fig_fontsize=12,suptit_=None,method='braycurtis'):
    
    
    """
    Plots several PCoA plots in one figure
    
    Parameters
    ----------
    
    otupcoadfs : list of pandas dataframes
        Dataframes to plot
    
    mapping : pandas dataframe
        A mapping dataframe containing the values to plot
    
    catvis : str
        The column name in the mapping file to plot
    
    fig_fontsize : int
        font size of the figure
    
    suptit_ : (False) Bool or str
        If str is passed then there will be a sup title
    
    method : int
        dissimilarity metric
    
    
    Returns
    -------
    matplotlib figure
    
    """
    
    if  all(isinstance(item, str) for item in list(mapping.T[catvis])) or all(isinstance(item, bool) for item in list(mapping.T[catvis])):
        if len(otupcoadfs)==1:
            fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(8, 8)) 
            otupcoadf=otupcoadfs[0]
            title='Mean Sequences Per Sample:'+str(np.mean(otupcoadf.sum())).split('.')[0]
            pcaplot=pcoa(DistanceMatrix(pdist(otupcoadf.as_matrix().T,method),list(otupcoadf.columns))).samples[['PC1','PC2','PC3']]
            pcaplot[catvis]=list(mapping.T[catvis])
            for ((key, grp)) in pcaplot.groupby(catvis):
                ax1.scatter(grp['PC1'], grp['PC2'], color=next(ax1._get_lines.prop_cycler)['color'], label=key, s=50)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
            ax1.set_title(title,fontsize=fig_fontsize,y=1.05)
            ax1.set_ylabel('$PC-1$',fontsize=fig_fontsize-2)
            ax1.set_xlabel('$PC-2$',fontsize=fig_fontsize-2) 
            if suptit_!=None:
                plt.suptitle(suptit_,y=1.05,fontsize=fig_fontsize+8)
            return fig
            
        else:
            fig, axn = plt.subplots(ncols=len(otupcoadfs), nrows=1, figsize=(15, 5)) 

            for ax1,(count_,otupcoadf) in zip(axn.flat,enumerate(otupcoadfs)):
                title=str(np.mean(otupcoadf.sum())).split('.')[0]+' Sequences Per Sample'
                pcaplot=pcoa(DistanceMatrix(pdist(otupcoadf.as_matrix().T,method),list(otupcoadf.columns))).samples[['PC1','PC2','PC3']]
                pcaplot[catvis]=list(mapping.T[catvis])
                for ((key, grp)) in pcaplot.groupby(catvis):
                    ax1.scatter(grp['PC1'], grp['PC2'], color=next(ax1._get_lines.prop_cycler)['color'], label=key, s=50)
                ax1.set_xticks([])
                ax1.set_yticks([])
                if count_==len(otupcoadfs)-1:
                    ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
                ax1.set_title(title,fontsize=fig_fontsize,y=1.05)
                ax1.set_ylabel('$PC-1$',fontsize=fig_fontsize-2)
                ax1.set_xlabel('$PC-2$',fontsize=fig_fontsize-2) 
            if suptit_!=None:
                plt.suptitle(suptit_,y=1.05,fontsize=fig_fontsize+8)
            return fig
    else:
        if len(otupcoadfs)==1:
            otupcoadf=otupcoadfs[0]
            pcaplot=pcoa(DistanceMatrix(pdist(otupcoadf.as_matrix().T,method),list(otupcoadf.columns))).samples[['PC1','PC2','PC3']]
            pcaplot[catvis]=list(mapping.T[catvis])
            fig, (ax1) = plt.subplots(ncols=1, nrows=1, figsize=(8, 6),sharey=False)
            pcaplot.plot.scatter(x='PC1', y='PC2', c=str(catvis),cmap='RdBu_r', s=50,ax=ax1);
            ax1.set_xlabel('$PC-1$')
            ax1.set_ylabel('$PC-2$')
            ax1.set_xticks([])
            ax1.set_yticks([])
            return fig 
        else:
            fig, axn = plt.subplots(ncols=len(otupcoadfs), nrows=1, figsize=(15, 5)) 
            for ax1,(count_,otupcoadf) in zip(axn.flat,enumerate(otupcoadfs)):
                title=str(np.mean(otupcoadf.sum())).split('.')[0]+' Sequences Per Sample'
                pcaplot=pcoa(DistanceMatrix(pdist(otupcoadf.as_matrix().T,method),list(otupcoadf.columns))).samples[['PC1','PC2','PC3']]
                pcaplot[catvis]=list(mapping.T[catvis])
                pcaplot.plot.scatter(x='PC1', y='PC2', c=str(catvis),cmap='RdBu_r', s=50,ax=ax1);
                ax1.set_xlabel('$PC-1$')
                ax1.set_ylabel('$PC-2$')
                ax1.set_xticks([])
                ax1.set_yticks([])
                ax1.set_title(title,fontsize=fig_fontsize,y=1.05)
            return fig 

def mean_KL(a,b):
    
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
    a=np.array(a.copy())
    b=np.array(b.copy())
    for i in range(a.shape[0]):
        kl += [kl_div(a[i],b[i])]
    return np.mean(kl)         

    

def build_block_model(rank,hoced,hsced,spar,C_,num_samples,num_features,overlap=0,mapping_on=True):
    
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
    X_true=block_diagonal_gaus(num_samples,num_features,rank,overlap,minval=.01,maxval=C_)
    if mapping_on==True:
    # make a mock mapping data
        mappning_=pd.DataFrame(np.array([['Cluster %s'%str(x)]*int(num_samples/rank) for x in range(1,rank+1)]).flatten()
                               ,columns=['example'],index=['sample_'+str(x) for x in range(0,num_samples-2)])

    X_noise=X_true.copy()
    X_noise=np.array(X_noise)
    # add Homoscedastic noise
    X_noise=Homoscedastic(X_noise,hoced)
    # add Heteroscedastic noise
    X_noise=Heteroscedastic(X_noise,hsced)
    # Induce low-density into the matrix
    X_noise=Subsample(X_noise,spar,num_samples)

    #return the base truth and noisy data
    if mapping_on==True:
        return X_true,X_noise,mappning_
    else:
        return X_true,X_noise

def minimize_model(x0,bnds,X_true_):
    
    def add_noise_min(x_o,X_true=X_true_):

        rank=3
        hoced=x_o[1]
        hsced=x_o[2]
        spar=x_o[3]
        C_=x_o[4]
        overlap_=x_o[5]

        num_samples=int(X_true.shape[0])
        num_features=int(X_true.shape[1])

        # make a mock OTU table
        X_true_=block_diagonal_gaus(num_samples,num_features,rank,overlap=overlap_,minval=.01,maxval=C_)
        # make a mock mapping data
        mappning_=pd.DataFrame(np.array([['Cluster %s'%str(x)]*int(num_samples/rank) for x in range(1,rank+1)]).flatten()
                               ,columns=['example'],index=['sample_'+str(x) for x in range(0,num_samples-2)])

        X_noise=X_true.copy()
        X_noise=np.array(X_noise)
        # add Homoscedastic noise
        X_noise=Homoscedastic(X_noise,hoced)
        # add Heteroscedastic noise
        X_noise=Heteroscedastic(X_noise,hsced)
        # Induce low-density into the matrix
        X_noise=Subsample(X_noise,spar,num_samples)

        #return the noisy data
        return mean_squared_error(X_true.T,X_noise) 
    

    model_fit=minimize(add_noise_min,x0, bounds=bnds)
    return model_fit


def build_grad_model(hoced,hsced,spar,sigma,C_,num_samples,num_features):
    
    
    rank = 2**3 - 1
    gradient = np.linspace(0, 10, num_samples)
    mu = np.linspace(0, 10, num_features)
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    table = np.vstack(xs).T
    table = pd.DataFrame(table)
    table.index = table.index.astype(np.str)
    table.columns = table.columns.astype(np.str)
    table=table*C_
    X_true=table.values.T.copy()
    
    X_noise=X_true.copy()
    X_noise=np.array(X_noise)
    # add Homoscedastic noise
    X_noise=Homoscedastic(X_noise,hoced)
    # add Heteroscedastic noise
    X_noise=Heteroscedastic(X_noise,hsced)
    # Induce low-density into the matrix
    X_noise=Subsample(X_noise,spar,num_samples)
        
    # make a mock mapping data
    mappning_=pd.DataFrame(np.array([x for x in range(0,num_samples)]),columns=['Gradient'],index=['sample_'+str(x) for x in range(0,num_samples)])
    mappning_=mappning_.apply(pd.to_numeric, errors='ignore')

    #return the base truth and noisy data 
    return X_true,X_noise,mappning_


def minimize_model_grad(x0,bnds,X_true_):


    def add_noise_min_grad(x_o,X_true=X_true_):

        #build model and minmize kl-div

        hoced=x_o[0]
        hsced=x_o[1]
        spar=x_o[2]
        sigma = x_o[3]
        C_ = x_o[4]
        
        num_samples=int(X_true.shape[0])
        num_features=int(X_true.shape[1])

        rank = 2**3 - 1
        gradient = np.linspace(0, 10, num_samples)
        mu = np.linspace(0, 10, num_features)
        xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
              for i in range(len(mu))]
        table = np.vstack(xs).T
        table = pd.DataFrame(table)
        table.index = table.index.astype(np.str)
        table.columns = table.columns.astype(np.str)
        table=table*C_

        X_noise=table.values.T.copy()
        X_noise=np.array(X_noise)
        # add Homoscedastic noise
        X_noise=Homoscedastic(X_noise,hoced)
        # add Heteroscedastic noise
        X_noise=Heteroscedastic(X_noise,hsced)
        # Induce low-density into the matrix
        X_noise=Subsample(X_noise,spar,num_samples)

        #return the noisy data 
        return mean_squared_error(X_true.T,X_noise) 
    
    model_fit=minimize(add_noise_min_grad,x0, bounds=bnds)
    
    return model_fit












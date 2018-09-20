import edward as ed
import tensorflow as tf

from gneiss.cluster import random_linkage
from edward.models import Normal, Poisson
from edward.models import PointMass

from biom import Table
import pandas as pd
import numpy as np

from numpy.linalg import matrix_rank
from scipy.sparse import coo_matrix
from skbio.stats.composition import _gram_schmidt_basis, ilr, clr_inv
from scipy.stats import norm

# helper function to perform sparse matrix multiplication
def sparse_matmul(A, B, row_index, col_index):
    """ Sparse matrix multiplication.
    
    This will try to evaluate the following product
    
    A[i] @ B[j]
    
    where i, j are the row and column indices specified in `indices`.
    
    Parameters
    ----------
    A : tf.Tensor
    Left 2D tensor
    B : tf.Tensor
    Right 2D tensor
    row_idx : tf.Tensor
    Row indexes to access in A
    col_idx : tf.Tensor
    Column indexes to access in B
    
    
    Returns
    -------
    tf.Tensor
    Result stored in a sparse tensor format, where the values
    are derived from A[i] @ B[j] where i, j are the row and column
    indices specified in `indices`.
    """
    A_flat = tf.gather(A, row_index, axis=0)
    B_flat = tf.transpose(tf.gather(B, col_index, axis=1))
    values = tf.reduce_sum(tf.multiply(A_flat, B_flat), axis=1)
    return values


def get_batch(M, Y):
    """ Get's batch data
    Parameters
    ----------
    M : int
    batch size
    Y : scipy.sparse.coo_matrix
    Scipy sparse matrix in COO-format.
    Returns
    -------
    batch_row : np.array
    Selected rows
    batch_col : np.array
    Selected columns
    batch_data : np.array
    Selected data
    """
    halfM = M // 2
    y_data = Y.data
    y_row = Y.row
    y_col = Y.col
    # get positive sample
    positive_idx = np.random.choice(len(y_data), halfM)
    positive_row = y_row[positive_idx]
    positive_col = y_col[positive_idx]
    positive_data = y_data[positive_idx]
    
    # store all of the positive (i, j) coords
    idx = np.vstack((y_row, y_col)).T
    idx = set(map(tuple, idx.tolist()))
    
    # get negative sample
    N, D = Y.shape
    negative_row = np.zeros(halfM)
    negative_col = np.zeros(halfM)
    negative_data = np.zeros(halfM)
    for k in range(halfM):
        i, j = np.random.randint(N), np.random.randint(D)
        while (i, j) in idx:
            i, j = np.random.randint(N), np.random.randint(D)
        negative_row[k] = i
        negative_col[k] = j
    batch_row = np.hstack((positive_row, negative_row))
    batch_col = np.hstack((positive_col, negative_col))
    batch_data = np.hstack((positive_data, negative_data))
    return batch_row, batch_col, batch_data

def ilr_to_clr(X, psi):
    """ Sparse ilr to clr transform.
    This is also known as the log_softmax transform with an extra
    rotation transform specified by the orthonormal matrix `psi`.
    Parameters
    ----------
    X : tf.Tensor
    Tensor of logits with dimensions N x D-1 to be transformed
    psi : tf.SparseTensor
    SparseTensor of dimensions D x D-1 used to transform `X`
    into a N x D dimesional tensor of logits.
    Returns
    -------
    tf.Tensor
    A N x D dimensional tensor of logits.  This also refered to as
    clr coordinates.
    """
    return tf.transpose(tf.sparse_tensor_dense_matmul(psi, tf.transpose(X)))


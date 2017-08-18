import numpy as np
from numpy import inf
from numpy.matlib import repmat
from numpy.linalg import norm
from scipy.sparse.linalg import svds
from skbio.stats.composition import ilr, clr_inv, _gram_schmidt_basis


def coptspace(M_E, r, niter, tol):
    """ Modified optspace algorithm for compositional data.

    Parameters
    ----------
    M_E, r, niter, tol

    Returns
    -------
    X, S, Y :
       SVD factors in ilr space

    """
    # ilr basis for performing the standard ilr transform.
    basis = _gram_schmidt_basis(M_E.shape[1])
    # convert everything to logs for ilr transform.
    # note - there will be -infs!
    M_E = np.log(M_E)
    # perform mean imputation scheme, so that we are only
    # measuring the subcompositional stress.
    M_E = np.vstack([P_E(M_E[i, :]) for i in range(M_E.shape[0])])
    # perform the actual ilr transform
    M_E = basis.dot(M_E)

    # TODO: will need to be careful here about setting things to zero
    E = np.abs(M_E) > 1e-10
    rescal_param = np.sqrt( np.count_nonzero(E) * r / norm(M_E, 'fro') ** 2 )
    M_E = M_E * rescal_param ;

    # TODO: Add in trimming for rows and columns
    X0, S0, Y0 = svds(M_E, r)
    n, m = M_E.shape
    nnz = np.count_nonzero(E)
    eps = nnz / np.sqrt(m*n)
    X0 = X0 * np.sqrt(n)
    Y0 = Y0 * np.sqrt(m)
    S0 = S0 / eps

    # wtf does 10000 come from?
    m0 = 10000

    rho = eps * n
    X, Y = X0, Y0.T;

    S = getoptS(X, Y, M_E, E)

    # TODO: Will need to figure out what norm is appropriate for termination.
    ft = M_E - X.dot(S).dot(Y.T)
    dist = np.zeros(niter + 1)
    dist[0] = norm( np.multiply(ft, E) ,'fro') / np.sqrt(nnz)
    for i in range(niter):
        W, Z = gradFt(X, Y, S, M_E, E, m0, rho);

        # Line search for the optimum jump length
        t = getoptT(X, W, Y, Z, S, M_E, E, m0, rho) ;
        X = X + t*W
        Y = Y + t*Z

        S = getoptS(X, Y, M_E, E) ;

        # Compute the distortion
        ft = M_E - X.dot(S).dot(Y.T)
        dist[i+1] = norm( np.multiply(ft, E) ,'fro') /  np.sqrt(nnz)
        if( dist[i+1] < tol ):
            break ;
    S = S /rescal_param ;

    return X, S, Y, dist


def optspace(M_E, r, niter, tol):
    """
    Parameters
    ----------
    M_E, r, niter, tol

    Returns
    -------
    X, S, Y
    """
    E = M_E > 0

    rescal_param = np.sqrt( np.count_nonzero(E) * r / norm(M_E, 'fro') ** 2 ) ;
    M_E = M_E * rescal_param ;
    # TODO: Add in trimming for rows and columns
    X0, S0, Y0 = svds(M_E, r)
    n, m = M_E.shape
    nnz = np.count_nonzero(E)
    eps = nnz / np.sqrt(m*n)
    X0 = X0 * np.sqrt(n)
    Y0 = Y0 * np.sqrt(m)
    S0 = S0 / eps

    # wtf does 10000 come from
    m0 = 10000

    rho = eps * n
    X, Y = X0, Y0.T;

    S = getoptS(X, Y, M_E, E)
    ft = M_E - X.dot(S).dot(Y.T)
    dist = np.zeros(niter + 1)
    dist[0] = norm( np.multiply(ft, E) ,'fro') / np.sqrt(nnz)
    for i in range(niter):
        W, Z = gradF_t(X, Y, S, M_E, E, m0, rho);

        # Line search for the optimum jump length
        t = getoptT(X, W, Y, Z, S, M_E, E, m0, rho) ;
        X = X + t*W
        Y = Y + t*Z

        S = getoptS(X, Y, M_E, E) ;

        # Compute the distortion
        ft = M_E - X.dot(S).dot(Y.T)
        dist[i+1] = norm( np.multiply(ft, E) ,'fro') /  np.sqrt(nnz)
        if( dist[i+1] < tol ):
            break ;
    S = S /rescal_param ;
    return X, S, Y, dist

def F_t(X, Y, S, M_E, E, m0, rho):
    """
    Parameters
    ----------
    X, Y, S, M_E, E, m0, rho

    Notes
    -----
    M ~ XSY
    """
    n, r = X.shape
    out1 = np.sum(
        np.sum((
            np.multiply((X.dot(S).dot(Y.T) - M_E), E)**2))) / 2
    out2 =  rho * G(Y, m0, r)
    out3 =  rho * G(X, m0, r)
    out = out1 + out2 + out3
    return out

def G(X, m0, r):
    """
    Parameters
    ----------
    X, m0, r
    """
    z = np.sum(X**2, axis=1) / (2 * m0 * r) ;
    y = np.exp( (z-1)**2 ) - 1 ;
    y[z<1] = 0
    y[y == np.inf] = 0
    return y.sum()

def gradF_t(X, Y, S, M_E, E, m0, rho):
    """
    Parameters
    ----------
    X, Y, S, M_E, E, m0, rho

    Returns
    -------
    W, Z
    """
    n, r = X.shape
    m, r = Y.shape

    XS = X.dot(S)
    YS = Y.dot(S.T)
    XSY = XS.dot(Y.T)

    Qx = X.T.dot( np.multiply((M_E - XSY), E)).dot(YS) / n;
    Qy = Y.T.dot( np.multiply((M_E - XSY), E).T).dot(XS) / m;
    W = np.multiply((XSY - M_E), E).dot(YS) + X.dot(Qx) + rho * Gp(X, m0, r);
         Z = np.multiply((XSY - M_E), E).T.dot(XS) + Y.dot(Qy) + rho * Gp(Y, m0, r);
    return W, Z

def Gp(X, m0, r):
    """
    X, m0, r
    """
    z = np.sum(X**2, axis=1) /(2 * m0 * r)
    z = 2 * np.multiply(np.exp( (z-1)**2 ), (z-1))

    z[z < 0] = 0
    z = z.reshape(len(z), 1)
    out = np.multiply(X, repmat(z, 1, r)) / (m0 * r)
    return out

def P_E(x):
    """ Replaces all infs with the running mean.

    Parameters
    ----------
    x : np.array
       Compositions with some zeros.

    Returns
    -------
    x_: np.array
       Imputed compositions
    """
    x_ = x.copy()
    for i in range(len(x)):
        if x[i] == -inf:
            x_[i] = x[:i].mean()
    return x_

def getoptT(X, W, Y, Z, S, M_E, E, m0, rho):
    """ X, W, Y, Z, S, M_E, E, m0, rho """
    norm2WZ = norm(W, 'fro')**2 + norm(Z, 'fro')**2

    # this is the resolution limit (t > 2**-20
    n_intervals = 20
    f = np.zeros(n_intervals+1)
    f[0] = F_t(X, Y, S, M_E, E, m0, rho)
    t = -1e-1

    for i in range(n_intervals):

        f[i+1] = F_t(X+t*W, Y+t*Z, S, M_E, E, m0, rho)
        if( (f[i+1] - f[1]) <= .5*t*norm2WZ ):
            return t
        t = t/2
    return t


def getoptS(X, Y, M_E, E):
    """ X, Y, M_E, E """
    n, r = X.shape

    C = np.ravel(X.T.dot(M_E).dot(Y))
    A = np.zeros((r*r, r*r))

    for i in range(r):
        for j in range(r):
            ind = j*r + i
            temp = np.multiply(X[ i,:].dot(Y[j,:].T), E)
            temp = X.T.dot(temp).dot(Y)
            A[:, ind] = np.ravel(temp)
    S = np.linalg.lstsq(A, C,
                        rcond = np.finfo(np.double).eps*int(max(A.shape)))[0]

    S = S.reshape((r, r))
    return S

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
from scipy import sparse
from DEICODE.ilr_pca import sparse_matmul
from DEICODE.ilr_pca import get_batch
from DEICODE.ilr_pca import ilr_to_clr

class ILR_PCA(object):
    
    
    """class for ILR-PCA on Edward and TensorFlow
    
    Requirements:
        TensorFlow version >= 0.6
        Edward version >= 1.3.5
        numpy>=1.7
        six>=1.10.0
        
    """
    
    def __init__(self, table, metadata,iterations=1000,learning_rate = 1e-1,batch_size = 1000
                 ,count_model='Poisson',alpha_mean = 0,alpha_scale = None
                 ,theta_mean = 0,theta_scale = None,gamma_mean = 0
                 ,gamma_scale = None,beta_mean = 0,beta_scale = 1):
        
        if batch_size>iterations:
            batch_size=iterations
        if iterations>batch_size:
            iterations=batch_size
        
        table=sparse.csr_matrix(np.array(table).copy())

        if alpha_scale is None:
            alpha_scale=2*np.sqrt(table.todense().mean()/matrix_rank(table.todense()))
        if theta_scale is None:
            theta_scale=2*np.sqrt(table.todense().mean()/matrix_rank(table.todense()))
        if gamma_scale is None:
            gamma_scale=2*np.sqrt(table.todense().mean()/matrix_rank(table.todense()))
        if beta_scale is None:
            beta_scale=2*np.sqrt(table.todense().mean()/matrix_rank(table.todense()))

        self.iterations=iterations
        self.table = table
        self.metadata = metadata
        self.batch_size=batch_size
        self.count_model=count_model
        self.learning_rate=learning_rate
        
        self.D, self.N = table.shape    
        self.p = metadata.shape[1]   # number of covariates
        
        self.basis = coo_matrix(_gram_schmidt_basis(self.D), dtype=np.float32).T

        # dummy variables for mini-batch size
        self.batch_row = tf.placeholder(tf.int32, shape=[batch_size], name='batch_i')
        self.batch_col = tf.placeholder(tf.int32, shape=[batch_size], name='batch_j')

        # global bias
        self.alpha = Normal(loc=tf.zeros([]) + alpha_mean,
                       scale=tf.ones([]) * alpha_scale,
                       name='alpha')
        # sample bias                                                                                    
        self.theta = Normal(loc=tf.zeros([self.N, 1]) + theta_mean,
                       scale=tf.ones([self.N, 1]) * theta_scale,
                       name='theta')
        # species bias
        self.gamma = Normal(loc=tf.zeros([1, self.D-1]) + gamma_mean,
                       scale=tf.ones([1, self.D-1]) * gamma_scale, 
                       name='gamma')

        # Specify sample principal axes
        self.U = Normal(loc=tf.zeros([self.N, self.p]),
                   scale=tf.ones([self.N, self.p]),
                   name='U')
        self.Uprime = tf.concat([self.theta, tf.ones([self.N, 1]), self.U], axis=1)

        # Specify feature principal axes
        self.V = Normal(loc=tf.zeros([self.p, self.D-1]),
                   scale=tf.ones([self.p, self.D-1]), 
                   name='V')
        self.Vprime = tf.concat([tf.ones([1, self.D-1]), self.gamma, self.V], axis=0)

        # Convert basis to SparseTensor
        self.psi = tf.SparseTensor(
            indices=np.mat([self.basis.row, self.basis.col]).transpose(),
            values=self.basis.data,
            dense_shape=self.basis.shape)

        # clr transform coefficients first                                                               
        self.cV = ilr_to_clr(self.Vprime, self.psi)

        # retrieve entries selected by index
        self.eta = sparse_matmul(self.Uprime, self.cV, 
            row_index=self.batch_row, col_index=self.batch_col)

        if count_model == "Poisson":
            self._build_Poisson_algorithm()
        else:
            raise ValueError("The attribute count_model must be in {'Poisson'} (Possble others in future versions)")
        
        return
    
    
    def _build_Poisson_algorithm(self):
        
        """build dataflow graph for Poisson """

        # obtain counts                                          
        self.Y = Poisson( rate=tf.exp(self.eta + self.alpha), name='Y' ) 

        # initialize the summaries
        summary_dir = 'summary-dir'

        with tf.variable_scope(tf.get_variable_scope(), reuse=tf.AUTO_REUSE):

            # These are the posterior distributions.
            tf.set_random_seed(0)

            qalpha_vars = [tf.get_variable("qalpha/loc", []),
                           tf.get_variable("qalpha/scale", [])]
            self.qalpha = Normal(loc=qalpha_vars[0],
                            scale=tf.nn.softplus(qalpha_vars[1]))

            qgamma_vars = [tf.get_variable("qgamma/loc", [1, self.D-1]),
                           tf.get_variable("qgamma/scale", [1, self.D-1])]
            self.qgamma = Normal(loc=qgamma_vars[0],
                            scale=tf.nn.softplus(qgamma_vars[1]))

            qtheta_vars = [tf.get_variable("qtheta/loc", [self.N, 1]),
                           tf.get_variable("qtheta/scale", [self.N, 1])]
            self.qtheta = Normal(loc=qtheta_vars[0],
                            scale=tf.nn.softplus(qtheta_vars[1]))

            qU_vars = [tf.get_variable("qU/loc", [self.N, self.p]),
                       tf.get_variable("qU/scale", [self.N, self.p])]
            self.qU = Normal(loc=qU_vars[0],
                        scale=tf.nn.softplus(qU_vars[1]))

            qV_vars = [tf.get_variable("qV/loc", [self.p, self.D-1]),
                       tf.get_variable("qV/scale", [self.p, self.D-1])]
            self.qV = Normal(loc=qV_vars[0],
                        scale=tf.nn.softplus(qV_vars[1]))

            # a placeholder for the microbial counts
            # since we will be manually feeding it into the inference via minibatch SGD
            self.Y_ph = tf.placeholder(tf.float32, shape=[self.batch_size], name='Y_placeholder')

            self.inference_u = ed.KLqp({
                self.U: self.qU,
                self.theta: self.qtheta,
                self.alpha: self.qalpha,
                self.gamma: self.qgamma}, 
                data={self.Y: self.Y_ph, self.V: self.qV})
            self.inference_v = ed.KLqp({
                self.V: self.qV,
                self.theta: self.qtheta,
                self.alpha: self.qalpha,
                self.gamma: self.qgamma}, 
                data={self.Y: self.Y_ph, self.U: self.qU})
            
            
            self.inference_u.initialize(var_list=qU_vars + qalpha_vars + qtheta_vars,
                                   n_samples=5,                       
                                   logdir=summary_dir)

            self.inference_v.initialize(var_list=qV_vars + qalpha_vars + qgamma_vars,
                                   n_samples=5,
                                   logdir=summary_dir, optimizer="adadelta")

            # adds checks for nans
            #tf.add_check_numerics_ops() so slow!
            
            return
        
        
    def complete(self,minval=None,maxval=None,clip_high=True):
        
        if maxval=='mean_iter':
            x_=self.table.todense()
            maxval=x_.T.sum(1)/(x_.T!=0).sum(1)
        if maxval=='mean':
            x_=self.table.todense()
            maxval=np.mean(x_.T.sum(1)/(x_.T!=0).sum(1))

        # initialize all tensorflow variables
        sess = ed.get_session()
        tf.global_variables_initializer().run()

        if self.count_model == "Poisson":
            return self._run_Poisson(sess,minval,maxval,clip_high)
        else:
            raise ValueError
        return

    def _run_Poisson(self, sess,minval,maxval,clip_high):
        
        losses = np.array([0.] * self.iterations)
        errors = np.array([0.] * self.iterations)
        data = self.table.tocoo().T
        y_row, y_col, y_data = data.row, data.col, data.data
        
        #from bayesian_regression.util.inference import get_batch
        for i in range(self.iterations):
            if not np.any(np.isnan(sess.run(self.qU.mean()))):
                u_ = sess.run(self.qU.mean())
                v_ = sess.run(self.qV.mean())
                theta_ = sess.run(self.qtheta.mean())
                gamma_ = sess.run(self.qgamma.mean())
                alpha_ = sess.run(self.qalpha.mean())

                log_u = ((u_ @ v_ + theta_).T + gamma_.T + alpha_).T @ self.basis.T
                err = np.mean((np.exp(log_u[y_row, y_col]) - y_data).ravel() ** 2)
                errors[i] = err

            # get batches
            idx_row, idx_col, idx_data = get_batch(M=self.batch_size, Y=data)

            self.inference_u.update(
                feed_dict={self.batch_row: idx_row, self.batch_col: idx_col, self.Y_ph: idx_data})
            info_dict = self.inference_v.update(
                feed_dict={self.batch_row: idx_row, self.batch_col: idx_col, self.Y_ph: idx_data})
            self.inference_v.print_progress(info_dict)
            losses[i] = info_dict['loss']

        u_ = sess.run(self.qU.mean())
        v_ = sess.run(self.qV.mean())
        theta_ = sess.run(self.qtheta.mean())
        gamma_ = sess.run(self.qgamma.mean())
        alpha_ = sess.run(self.qalpha.mean())

        log_u = ((u_ @ v_ + theta_).T + gamma_.T + alpha_).T @ self.basis.T

        etable = np.zeros((self.N, self.D))
        for i in range(etable.shape[0]): 
            etable[i] = np.random.poisson(np.exp(log_u[i])) #completed table
        etable=etable.T
        
        otu_org=self.table.todense()
        #clip high values
        if clip_high==True:
            shape=list(etable.shape)
            for x in range(0, shape[0]):
                for y in range(0, shape[1]):
                    #clip values where imoute is higher than orig
                    if otu_org[x, y]!=0 and etable[x, y]>otu_org[x, y]:
                        etable[x, y]=otu_org[x, y]
                    #after impute that were there before (should not happen)
                    if otu_org[x, y]!=0 and etable[x, y]==0:
                        etable[x, y]=otu_org[x, y]

        #clip high values with iter
        if maxval is not None and isinstance(maxval, (list, tuple, np.ndarray)):
            shape=list(etable.shape)
            for x in range(0, shape[0]):
                for y,maxval_ in zip(range(0, shape[1]),maxval):
                    if otu_org[x, y]==0 and etable[x, y]>maxval_:
                        etable[x, y]=maxval_

        #clip high values with float (if both None skip)
        if maxval is not None and isinstance(maxval, (int, float)):
            shape=list(etable.shape)
            for x in range(0, shape[0]):
                for y in range(0, shape[1]):
                    if otu_org[x, y]==0 and etable[x, y]>maxval:
                        etable[x, y]=maxval

        #clip values below min
        if minval!=None:
            etable[etable<minval]=minval

        return etable


import tensorflow as tf
import edward as ed
from edward.models import Poisson, Gamma, Empirical, PointMass
from numpy.linalg import matrix_rank
import numpy as np
import pandas as pd

class GaPNMF(object):
    
        
    """class for Gamma-Poisson factorization GaP NNMF on Edward and TensorFlow
    
    Requirements:
        TensorFlow version >= 0.6
        Edward version >= 1.3.5
        numpy>=1.7
        six>=1.10.0
        
    """

    def __init__(self,X,rank=None,n_iter=500,t=500,alpha=None,beta=None,gamma=None,count_model="Gamma-Poisson"):

        #Empirical approximation cannot hold any more samples, so t must be > n_iter
        if n_iter>t:
            n_iter=t
            
        if rank==None:
            rank=int(matrix_rank(X))
        else:
            rank=rank

        # get missing data
        self.count_model=count_model
        self.n_iter=n_iter
        self.X=np.array(X).copy()
        
        idx = X.nonzero()
        tidx = tf.constant(np.column_stack(idx))
        y = X[idx]
        n,m = X.shape
        y = np.ceil(y)

        if self.count_model == "Gamma-Poisson":
            
            #priors
            # concentration (alpha) , rate (beta)
            mean=np.mean(y)
            stddev=np.std(y)
            if alpha==None:
                alpha = np.float32((mean / stddev)**2) #alpha
            if beta==None:
                beta = np.float32(mean / stddev**2)  #beta
            if gamma==None:
                gamma=np.float32(beta/10) #pop betas
            
            act = Gamma(alpha, beta, sample_shape=n) # Sample activity
            pref = Gamma(alpha, act, sample_shape=rank) # Sample preference
            pop = Gamma(gamma, gamma, sample_shape=m) # OTU popularity
            attr = Gamma(alpha, pop, sample_shape=rank) # OTU attribute
            like = Poisson(tf.gather_nd(tf.matmul(pref, attr, transpose_a=True), tidx))

            #prost
            qact = Empirical(
                tf.nn.softplus(tf.Variable(tf.random_normal([t,n]))),)
            qpref = PointMass(
                tf.nn.softplus(tf.Variable(tf.random_normal([rank,n]))),)
            qpop = Empirical(
                tf.nn.softplus(tf.Variable(tf.random_normal([t,m]))),)
            qattr = PointMass(
                tf.nn.softplus(tf.Variable(tf.random_normal([rank,m]))),)

            #inference
            self.inference_e = ed.Gibbs(
                {act:qact, pop:qpop}, 
                data={like:y, pref:qpref, attr:qattr},)

            self.inference_m = ed.MAP(
                {pref:qpref, attr:qattr},
                data={like:y, act:qact, pop:qpop},)
            
            self.qpref=qpref
            self.qattr=qattr
            
            
            self._build_GaP_algorithm()
            
        else:
            raise ValueError("The attribute count_model must be in {'Poisson'} (Possble others in future versions)")

        return


    def _build_GaP_algorithm(self):

        self.inference_e.initialize()
        self.inference_m.initialize(n_iter=self.n_iter, optimizer="adadelta")

        return 

    def complete(self,minval=None,maxval='mean',clip_high=True):
        
        if maxval=='mean_iter':
            maxval=self.X.T.sum(1)/(self.X.T!=0).sum(1)
        if maxval=='mean':
            maxval=np.mean(self.X.T.sum(1)/(self.X.T!=0).sum(1))
        
        # initialize all tensorflow variables
        sess = ed.get_session()
        tf.global_variables_initializer().run()

        if self.count_model == "Gamma-Poisson":
            return self._run_Poisson(sess,minval,maxval,clip_high)
        else:
            raise ValueError
        return

    def _run_Poisson(self,sess,minval,maxval,clip_high):

        #pass in#
        loss=np.empty(self.n_iter,dtype=np.float32)
        for i in range(self.n_iter):
            info_dict_e=self.inference_e.update()
            info_dict_m=self.inference_m.update()
            loss[i]=info_dict_m["loss"]
            self.inference_m.print_progress(info_dict_m)

        pref = sess.run(self.qpref) # Infered sample preference.
        attr = sess.run(self.qattr) # Infered specia attribute.
        poisson = np.random.poisson

        etable=poisson(np.dot(pref.T, attr))

        otu_org=self.X
        
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


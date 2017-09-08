from __future__ import absolute_import, division, print_function
import unittest
import pandas as pd
import numpy as np
from gneiss import util
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
import numpy.testing as npt
import pandas.util.testing as pdt
from DEICODE.untangle import (machine_learning,biplot,relative_abund_,
                              encode_mapping,complete_matrix,features_ml)


class TestUntangle(unittest.TestCase):
    
    def test_biplot(self):
        
        
        n=10
        b=10
        r=2
        X_true=pd.DataFrame(util.block_diagonal(n,b,r),columns=[x for x in range(1,n+1)],index=[x for x in range(1,(b+1))])
        X_shuffle=X_true.T.sample(frac=1)
        X_shuffle=X_shuffle.T.sample(frac=1)
        bione,bitwo=biplot(X_shuffle,r=2)
        
        sample_groups_truth=[X_true.columns.tolist()[0:5] ,X_true.columns.tolist()[5:11]]
        OTU_groups_truth=[X_true.index.tolist()[0:5] ,X_true.index.tolist()[5:11]]
        
        OTU_groups=[list(np.sort([g_ for g_,key in zip(bione.index.tolist(),list(map(bool,bitwo.index.tolist()))) if not key])),
                    list(np.sort([g_ for g_,key in zip(bione.index.tolist(),list(map(bool,bitwo.index.tolist()))) if key]))]
            
            
        sample_groups=[list(np.sort([g_ for g_,key in zip(bione.columns.tolist(),list(map(bool,bitwo.columns.tolist()))) if key])),
                       list(np.sort([g_ for g_,key in zip(bione.columns.tolist(),list(map(bool,bitwo.columns.tolist()))) if not key]))]
        
        try:
            assert_array_equal(OTU_groups_truth,OTU_groups)
        except:
            assert_array_equal(OTU_groups_truth,[OTU_groups[1],OTU_groups[0]])
                
        try:
            assert_array_equal(sample_groups_truth,sample_groups)
        except:
            assert_array_equal(sample_groups_truth,[sample_groups[1],sample_groups[0]])


    def test_relative_abund_(self):
        
        truth=pd.DataFrame([[100.0,50.0],[0.0,50.0]],index=['o1','o2'],columns=['s1','s2'])
        test=relative_abund_(pd.DataFrame([[45,10],[0,10]],index=['o1','o2'],columns=['s1','s2']))
        assert_array_equal(truth,test)

    def test_encode_mapping(self):
        
        truth=pd.DataFrame([[0,0],[1,1]],index=['s1','s2'],columns=['env_1','env_2'])
        test=encode_mapping(pd.DataFrame([['A','C'],['B','D']],index=['s1','s2'],columns=['env_1','env_2']))
        assert_array_equal(truth,test)

    def test_complete_matrix(self):
        
        truth=np.array([[ 0.2,0.4,0.4,0.],[ 0.23683605,0.5,0.5,0.]])
        test=complete_matrix(np.array([[.2,.4,.4, 0],[0,.5,.5,0]]),rank=1,iteration=40)
        assert_array_almost_equal(truth,test,decimal=1)

    def test_features_ml(self):
        
        
        truth=pd.DataFrame([0.3,0.2,0],index=['o1','o2','o3'])
        truth_table=pd.DataFrame([[1000,0],[0,1000],[1000,1000]],index=['o1','o2','o3'],columns=['s1','s2'])
        truth_mapping=pd.DataFrame(['state_1','state_2'],index=['s1','s2'],columns=['state'])
        order_complete=features_ml(truth_table,truth_mapping,'state')['Importance']
        order_no_complete=features_ml(truth_table,truth_mapping,'state',complete=False)['Importance']
        assert_array_almost_equal(truth,pd.DataFrame(order_complete),decimal=1)
        assert_array_almost_equal(truth,pd.DataFrame(order_no_complete),decimal=1)

    def test_machine_learning(self):
        
        scores_truth=pd.DataFrame(np.array([[ 0.28,0.44899889],[ 0.72,0.44899889]]),index=['state_m1 (n=3) (labels=2)', 'state_m2 (n=3) (labels=2)'],columns=['Mean Matrix Completion (RF)', 'std Matrix Completion (RF)'])

        table=pd.DataFrame([[1000,0,0],[0,1000,0],[1000,0,0]],index=['o1','o2','o3'],columns=['s1','s2','s3'])
        mapping=pd.DataFrame(np.array([['state_1','state_2','state_1'],['state_1','state_2','state_2']]).T,index=['s1','s2','s3'],columns=['state_m1','state_m2'])
        scores_test,test_comp=machine_learning(table,mapping,mean_count=50)
        try:
            assert_array_almost_equal(np.array([[0.4,0.5],[0.6,0.5]]),scores_test,decimal=1)
        except:
            assert_array_almost_equal(np.array([[0.6,0.5],[0.4,0.5]]),scores_test,decimal=1)


if __name__ == '__main__':
    unittest.main()


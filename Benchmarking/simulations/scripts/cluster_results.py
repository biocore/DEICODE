#utils
import pandas as pd
import numpy as np
from gneiss.util import match
#KL-div
from scipy.stats import entropy
#PCoA and PERMANOVA
from scipy.spatial import distance
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
#Classifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
#RCPCA
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
from skbio.stats.composition import closure, clr, clr_inv

#import base truth data
base_truth=pd.read_csv('cluster_models/simulation_base_truth.csv', index_col=[0,1,2,3])
#import noisy and sparse data
subsampled=pd.read_csv('cluster_models/simulation_subsampled_noisy.csv', index_col=[0,1,2,3])
#intial dict to save results
results={}
ranks_ = set(subsampled.index.get_level_values('rank'))
#powers_ = set(subsampled.index.get_level_values('overlap'))
powers_ = [20]
depths_ = set(subsampled.index.get_level_values('sequence_depth'))
for count_,rank_ in enumerate(ranks_): #iter by rank (only two)
    for power_ in powers_: #iter by overlap between clusters
        for depth_ in depths_: #inter by seq. depth
            
            #get the data for that subset
            subtmp=subsampled.loc[(rank_,power_,depth_,),:].copy().T 
            #get the base truth data for that subset
            basetmp_sub=base_truth.loc[(rank_,power_,depth_,),:].copy().T
            # sub sampled
            subtmp_sub = subtmp.copy()
            #meta on cluster
            meta = np.array([1]*int(subtmp.shape[0]/2)+[2]*int(subtmp.shape[0]/2)).T
            meta = pd.DataFrame(meta,index=subtmp.index,columns=['group'])

            # test KL with rcl
            X_sparse=rclr().fit_transform(subtmp_sub.copy())
            U,s,V=OptSpace(iteration=1000).fit_transform(X_sparse)
            clr_res = clr_inv(np.dot(np.dot(U,s),V.T))
            # use just kl_div here because already closed 
            kl_clr = entropy(closure(basetmp_sub).T,clr_res.T).mean()
            results[(rank_,power_,depth_,'rclr','KL-Div')]=[kl_clr]
            
            # test KL without rclr
            X_spn = np.array(subtmp_sub.copy()).astype(float)
            X_spn[X_spn==0] = np.nan
            U_,s_,V_=OptSpace(iteration=1000).fit_transform(X_spn)
            res_raw = np.dot(np.dot(U_,s_),V_.T)
            res_raw[res_raw<=0]=1
            kl_raw = entropy(closure(basetmp_sub).T,closure(res_raw).T).mean()
            results[(rank_,power_,depth_,'Raw Counts','KL-Div')]=[kl_raw]
            
            # f-stat
            resfclr=permanova(DistanceMatrix(distance.cdist(U,U)),meta['group'])['test statistic']
            rawfres=permanova(DistanceMatrix(distance.cdist(U_,U_)),meta['group'])['test statistic']
            results[(rank_,power_,depth_,'rclr','F-Statistic')]=[resfclr]
            results[(rank_,power_,depth_,'Raw Counts','F-Statistic')]=[rawfres]
            
            # KNN
            for U_tmp,method in zip([U,U_],['rclr','Raw Counts']):
                pcoa_tmp=pcoa(DistanceMatrix(distance.cdist(U_tmp,U_tmp))).samples
                pcoa_tmp.index=subtmp_sub.index
                # split
                X_train, X_test, y_train, y_test = train_test_split(pcoa_tmp, meta['group'].ravel(), 
                                                                    test_size=.4, 
                                                                    stratify=meta['group'].ravel(), 
                                                                    random_state=42)
                knn = KNeighborsClassifier(n_neighbors=13).fit(X_train, y_train)
                acc_ = accuracy_score(y_test,knn.predict(X_test).astype(int))
                results[(rank_,power_,depth_,method,'KNN_Accuracy')]=[acc_]
            
#convert the results to df to make it easier to read and plot       
resdf = pd.DataFrame(results).T
resdf.index.names = ['Rank', 'Overlap','Sequence_Depth','Method','Metric']
resdf.columns=['value']
resdf = resdf.reset_index()
resdf.to_csv('cluster_models/results.csv')
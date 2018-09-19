import pandas as pd
import numpy as np

from sklearn.metrics import mean_squared_error
from DEICODE.utils import mean_KL

from DEICODE.utils import rclr
from DEICODE.opt_space import optspace
from fancyimpute import KNN, SoftImpute, IterativeImputer, BiScaler
from skbio.stats.composition import clr_inv
from skbio.stats.composition import closure
from skbio.stats.composition import clr

#import base truth data
base_truth=pd.read_csv('cluster_models/simulation_base_truth.csv', index_col=[0,1,2,3])
#import noisy and sparse data
subsampled=pd.read_csv('cluster_models/simulation_subsampled_noisy.csv', index_col=[0,1,2,3])

#intial dict to save results
results={}
for rank_ in set(subsampled.index.get_level_values('rank')): #iter by rank (only two)
    for power_ in set(subsampled.index.get_level_values('overlap')): #iter by overlap between clusters
        for depth_ in set(subsampled.index.get_level_values('sequence_depth')): #inter by seq. depth

            subtmp=subsampled.loc[(rank_,power_,depth_,),:].copy() #get the data for that subset
            subtmp_rclr=rclr(subtmp) # take the robust clr transform
            basetmp=base_truth.loc[(rank_,power_,depth_,),:].copy() #get the base truth data for that subset
            
            U, s, V, _ =optspace(subtmp_rclr,r=rank_, niter=5, tol=1e-5) #fit opt-space to the r-clr data
            optcomp=U.dot(s).dot(V.T) # recompose the completed matrix 
            
            X_noisefancy=subtmp.astype(float).as_matrix().copy() #save as float 
            X_noisefancy[X_noisefancy==0]=np.nan #convert zeros to nan for fancy impute             
            X_filled_knn = KNN(verbose=False,min_value=1).fit_transform(X_noisefancy) # KNN fill 
            X_filled_softimpute = SoftImpute(verbose=False,min_value=1).fit_transform(X_noisefancy) #SoftImpute
            X_filled_iter=IterativeImputer(verbose=False,min_value=1).fit_transform(X_noisefancy) #Iter SVD
            
            #calculate the KL-divergence and then fill the dict of results 
            results[(rank_,power_,depth_,'OptSpace','KL-Div')]=[mean_KL((optcomp),rclr(basetmp))]
            results[(rank_,power_,depth_,'Psuedo','KL-Div')]=[mean_KL(rclr(subtmp+1),rclr(basetmp))]
            results[(rank_,power_,depth_,'KNN','KL-Div')]=[mean_KL(rclr(X_filled_knn),rclr(basetmp))]
            results[(rank_,power_,depth_,'Soft','KL-Div')]=[mean_KL(rclr(X_filled_softimpute),rclr(basetmp))]
            results[(rank_,power_,depth_,'Iter','KL-Div')]=[mean_KL(rclr(X_filled_iter),rclr(basetmp))]
            
#convert the results to df to make it easier to read and plot       
results=pd.DataFrame(results).T
results.index.names = ['Rank', 'Overlap','Sequence_Depth','Method','Metric']
results.columns=['value']
results.to_csv('cluster_models/results.csv')
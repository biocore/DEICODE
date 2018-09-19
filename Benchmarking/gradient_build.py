#Import 
import warnings; warnings.simplefilter('ignore') #for PCoA warnings
import pandas as pd
import numpy as np

#import data
from biom import load_table
from skbio.stats import subsample_counts
#MOCK data generation
from gneiss.util import match
from gneiss.sort import niche_sort
from DEICODE.utils import block_diagonal_gaus
from DEICODE.utils import build_grad_model
from DEICODE.utils import minimize_model
from DEICODE.utils import minimize_model_grad
#compostional transform
from skbio.stats.composition import clr

# import observation data
in_biom='gradient_models/88soils.biom' #import biom file
table = load_table(in_biom)
read_filter = lambda val, id_, md: sum(val) > 1000
table.filter(read_filter, axis='sample')
otutabledf=table.to_dataframe()
otutabledf=otutabledf.T

# Get OTU to taxa match
taxonomy=table.metadata_to_dataframe('observation')
taxonomy.columns=['kingdom', 'phylum', 'class', 'order', 
                             'family', 'genus', 'species']
taxonomy['taxonomy'] = taxonomy[taxonomy.columns].apply(lambda x: ';'.join(x), axis=1)

#mapping import 
map_file='gradient_models/88soils.txt' #import metadata
mappingdf= pd.read_table('%s'%map_file, index_col=0,low_memory=False)
mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)
mappingdf.index=list(map(str,mappingdf.index))
mappingdf=mappingdf.astype(str)

#match the tables
otutabledf,mappingdf=match(otutabledf,mappingdf)

otutabledf.columns=[str(x) for x in otutabledf.columns]
mappingdf=mappingdf.apply(pd.to_numeric, errors='ignore')

observed_table = niche_sort(otutabledf, mappingdf.ph)
mappingdf=mappingdf.T[observed_table.index].T
otutabledf=observed_table.copy()

otutabledf.to_dense().to_csv("gradient_models/base_model_soils_table.csv",sep=',', encoding='utf-8')
mappingdf.to_dense().to_csv("gradient_models/base_model_soils_meta.csv",sep=',', encoding='utf-8')


# hoced,hsced,spar,sigma,C_,num_samples,num_features
X_true=np.array(otutabledf.as_matrix())
x0 = [.5, .5, 1e2, 2, 1e2]
bnds = ((0,1),(0,1),(1,1e2),(1,3),(1,1e2))
model_fit=minimize_model_grad(x0,bnds,X_true)

save_base=[]
save_sub=[]
for sigma_ in [1.0,1.3,1.7,2.0]:
    
    #sequencing depth
    seq_depth={100:1e2,1000:1e3,2000:2e3,5000:5e3,10000:1e4}    
    
    # sub sample 
    for sub_depth,sub_ in seq_depth.items():

        #run model with fit variables and new variants
        base_truth,X_noise_sparse,mapping=build_grad_model(model_fit.x[0],model_fit.x[1], sub_depth
                                                           , sigma_
                                                           , sub_depth
                                                    ,otutabledf.shape[0]
                                                      ,1000)
        base_truth=pd.DataFrame(base_truth
                                ,index=[(sigma_,sub_,'OTU_'+str(x)) for x in range(base_truth.shape[0])]
                               ,columns=['sample_'+str(x) for x in range(base_truth.shape[1])])
        
        X_noise_sub=pd.DataFrame(X_noise_sparse
                        ,index=[(sigma_,sub_,'OTU_'+str(x)) for x in range(X_noise_sparse.shape[0])]
                       ,columns=['sample_'+str(x) for x in range(X_noise_sparse.shape[1])])

        #for X_noise_subsampled in Subsamples_noisy:
        save_base.append(base_truth)
        save_sub.append(X_noise_sub)

for df_,loc_ in zip([save_base,save_sub]
                ,['simulation_base_truth','simulation_subsampled_noisy']):

    df_=pd.concat(df_,axis=0)
    df_.index=pd.MultiIndex.from_tuples(df_.index)
    df_.index.names = ['band_width','sequence_depth','OTUs']
    df_.to_csv('gradient_models/'+loc_+'.csv') #save both and finish
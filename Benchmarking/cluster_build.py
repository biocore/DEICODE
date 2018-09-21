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
from deicode.utils import block_diagonal_gaus
from deicode.utils import build_block_model
from deicode.utils import minimize_model
#compostional transform
from skbio.stats.composition import clr


# import observation data
in_biom='cluster_models/keyboard.biom' #import biom file
table = load_table(in_biom)
read_filter_s = lambda val, id_, md: sum(val) > 0
read_filter_f = lambda val, id_, md: sum(val) > 0
table=table.filter(read_filter_s, axis='sample')
table=table.filter(read_filter_f, axis='observation')
otutabledf=table.to_dataframe()
otutabledf=otutabledf.T
otutabledf.drop_duplicates(inplace=True)

# Get OTU to taxa match
taxonomy=table.metadata_to_dataframe('observation')
taxonomy.columns=['kingdom', 'phylum', 'class', 'order', 
                             'family', 'genus', 'species']
taxonomy['taxonomy'] = taxonomy[taxonomy.columns].apply(lambda x: ';'.join(x), axis=1)

#mapping import 
map_file='cluster_models/keyboard.txt' #import metadata
mappingdf= pd.read_table('%s'%map_file, index_col=0,low_memory=False)
mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)
mappingdf.index=list(map(str,mappingdf.index))
mappingdf=mappingdf.astype(str)
mappingdf=mappingdf[~mappingdf.index.duplicated(keep='first')]

#match the tables
otutabledf,mappingdf=match(otutabledf,mappingdf[mappingdf['host_subject_id'].isin(['M2','M3','M9'])])

otutabledf=otutabledf.T[otutabledf.sum()>0].T
otutabledf=otutabledf[otutabledf.T.sum()>0]
otutabledf.columns=[str(x) for x in otutabledf.columns]

sorting_map={'M9':2,'M2':3,'M3':1}

mappingdf['host_num']=[int(sorting_map[x]) for x in mappingdf['host_subject_id']]
mappingdf=mappingdf.apply(pd.to_numeric, errors='ignore')

#sort by niche 
observed_table = niche_sort(otutabledf, mappingdf['host_num'])
mappingdf=mappingdf.T[observed_table.index].T
otutabledf=observed_table.copy()

otutabledf.to_dense().to_csv("cluster_models/base_model_keyboard_table.csv",sep=',', encoding='utf-8')
mappingdf.to_dense().to_csv("cluster_models/base_model_keyboard_meta.csv",sep=',', encoding='utf-8')

######### build the model #########
x0 = [3, 20, 20, 1e2, 1e2,1e1]
bnds = ((3,3),(0,1e2),(0,2e3),(0,1e10),(0,5e1),(1,10))
model_fit=minimize_model(x0,bnds,np.array(otutabledf.T[:104].T.as_matrix()))

save_base=[]
save_sub=[]
for rank_,overlap_ in zip([2,2,2,2],
                             [0,5,10,20]):
    
    #subsample_points=np.logspace(2,4,4)
    seq_depth={100:6e2,1000:1e3,2000:2e3,5000:5e3,10000:1e4}
    for sub_,model_peram in seq_depth.items():
        
        #run model with fit variables and new variants
        base_truth,X_noise_sub=build_block_model(rank_,  model_fit.x[1], model_fit.x[2], model_peram, model_peram
                                                 ,200,4000,overlap=overlap_
                                                 ,mapping_on=False)
        base_truth=pd.DataFrame(base_truth
                                ,index=[(rank_,overlap_,sub_,'OTU_'+str(x)) for x in range(base_truth.shape[0])]
                               ,columns=['sample_'+str(x) for x in range(base_truth.shape[1])])
        
        X_noise_sub=pd.DataFrame(X_noise_sub
                        ,index=[(rank_,overlap_,sub_,'OTU_'+str(x)) for x in range(X_noise_sub.shape[0])]
                       ,columns=['sample_'+str(x) for x in range(X_noise_sub.shape[1])])

        #for X_noise_subsampled in Subsamples_noisy:
        save_base.append(base_truth)
        save_sub.append(X_noise_sub)

for df_,loc_ in zip([save_base,save_sub]
                ,['simulation_base_truth','simulation_subsampled_noisy']):

    df_=pd.concat(df_,axis=0)
    df_.index=pd.MultiIndex.from_tuples(df_.index)
    df_.index.names = ['rank', 'overlap','sequence_depth','OTUs']
    df_.to_csv('cluster_models/'+loc_+'.csv') #save both and finish
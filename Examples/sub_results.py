import pandas as pd
import numpy as np
import os
from skbio.stats.ordination import OrdinationResults
from gneiss.util import match
from biom.util import biom_open
from biom import load_table
from collections import Counter
import numpy as np
import subprocess
import shutil 

def make_safe_dir(tmp_w):
    if not os.path.exists(tmp_w):
        os.makedirs(tmp_w)
    else:
        shutil.rmtree(tmp_w)           
        os.makedirs(tmp_w)

case_study={}
case_study['Sponges']={}
case_study['Sleep_Apnea']={}
case_study['AG']={}
case_study['Sponges']['factor']='health_status'
case_study['Sleep_Apnea']['factor']='exposure_type'
case_study['AG']['factor']='types_of_plants'
distances={}
meta={}
ordinations={}
# empty dict to save each case study in 
dir_path = os.getcwd()
for dataset_ in [dir_ for dir_ in os.listdir('data/') if dir_!='.DS_Store']:
    distances[dataset_]={}
    meta[dataset_]={}
    ordinations[dataset_]={}
    subpath_=os.path.join('sub_sample','biom_tables_'+dataset_)
    for sub_set in [dir_ for dir_ in os.listdir(subpath_) if dir_!='.DS_Store']:

        depths_=list(map(int,sub_set.split('_')))
        fold_=depths_[0]
        Nsamp_=depths_[1]

        # set emtpy 
        distances[dataset_][(fold_,Nsamp_)]={}
        meta[dataset_][(fold_,Nsamp_)]={}
        ordinations[dataset_][(fold_,Nsamp_)]={}
        #metadata for subset
        meta_ = pd.read_table(os.path.join(subpath_,sub_set,'metadata.tsv'), index_col=0)
        meta_.index=meta_.index.astype(str)

        # each distance 
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'GUniFrac_alpha_one_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        index_me=list(set(table_.index)&set(table_.columns))
        index_me=list(set(index_me)&set(meta_.index))
        #reindex meta
        meta[dataset_][(fold_,Nsamp_)]['metadata'] = meta_.reindex(index=index_me)

        #reindex all others the same
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['GUniFrac_Alpha_One'] = table_
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'GUniFrac_alpha_half_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['GUniFrac_Alpha_Half'] = table_
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'GUniFrac_alpha_zero_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['GUniFrac_Alpha_Zero'] = table_
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'Jaccard_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['Jaccard'] = table_
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'Bray_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['Bray_Curtis'] = table_
        table_ = pd.read_table(os.path.join(subpath_,sub_set,'Robust_Aitchison_Distance.tsv'), 
                               index_col=0,low_memory=False)
        table_.index=table_.index.astype(str)
        table_.columns=table_.columns.astype(str)
        table_ = table_.reindex(index=index_me,columns=index_me)
        distances[dataset_][(fold_,Nsamp_)]['Robust_Aitchison'] = table_

        # ordination type file
        in_ord = os.path.join(subpath_,sub_set,'RPCA_Ordination.txt')
        # get loadings from ordination files
        ordinations[dataset_][(fold_,Nsamp_)]['RPCA_Samples'] = OrdinationResults.read(in_ord).samples
        ordinations[dataset_][(fold_,Nsamp_)]['RPCA_Features'] = OrdinationResults.read(in_ord).features

from skbio import DistanceMatrix
from skbio.stats.distance import permanova

both_perm_res={}
perm_res={}
perm_res_tmp={}
for dataset_,subs in distances.items():
    perm_res[dataset_]={}
    perm_res_tmp[dataset_]={}
    
    # run on subset 
    arr_run=np.sort(list(set([j for i,j in distances[dataset_].keys()])))
    idx = np.round(np.linspace(1, len(arr_run) - 1, 6)).astype(int)
    arr_run=arr_run[idx]
    
    for (fold_,Nsamp_),methods_ in subs.items():
        # run a subset
        if Nsamp_ not in arr_run:
            continue 
        meta_=meta[dataset_][(fold_,Nsamp_)]['metadata']
        if len(meta_.index)<Nsamp_:
            continue 
        perm_res[dataset_][(fold_,Nsamp_)]={}
        perm_res_tmp[dataset_][(fold_,Nsamp_)]={}
        for method,dist_tmp in methods_.items():
            if method in ['Bray_Curtis','GUniFrac_Alpha_One', 'Robust_Aitchison']:
                perm_res[dataset_][(fold_,Nsamp_)][method]={}
                dist_tmp=DistanceMatrix(dist_tmp)
                perm_tmp=permanova(dist_tmp,meta[dataset_][(fold_,Nsamp_)]['metadata'][case_study[dataset_]['factor']].values)
                perm_res[dataset_][(fold_,Nsamp_)][method]['test statistic']=perm_tmp['test statistic']
                perm_res_tmp[dataset_][(fold_,Nsamp_)]=pd.DataFrame(perm_res[dataset_][(fold_,Nsamp_)])

    both_perm_res[dataset_]=pd.concat(perm_res_tmp[dataset_])

#PCoA
import warnings; warnings.simplefilter('ignore') #for PCoA warnings
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
#Classifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

both_nn={}
nn_res={}
nn_res_tmp={}
for dataset_,subs in distances.items():
    nn_res[dataset_]={}
    nn_res_tmp[dataset_]={}
    
    # run on subset 
    arr_run=np.sort(list(set([j for i,j in distances[dataset_].keys()])))
    idx = np.round(np.linspace(1, len(arr_run) - 1, 6)).astype(int)
    arr_run=arr_run[idx]

    # set splits for each dataset here
    max_size=np.max([j for i,j in distances[dataset_].keys()])
    dist_tmp=distances[dataset_][(1,max_size)]['Bray_Curtis'] #any method to split
    meta_=meta[dataset_][(1,max_size)]['metadata']
    pcoa_tmp=pcoa(DistanceMatrix(dist_tmp)).samples
    pcoa_tmp.index=dist_tmp.index
    # prepare binary classifier 
    y=meta_[case_study[dataset_]['factor']].values
    proc=preprocessing.LabelBinarizer()
    proc.fit(y)
    meta_['y_encode']=proc.transform(y)
    # split
    X_train, X_test, y_train, y_test = train_test_split(pcoa_tmp, meta_['y_encode'].ravel(), 
                                                        test_size=.4, 
                                                        stratify=meta_['y_encode'].ravel(), 
                                                        random_state=42)
    train_index =X_train.index
    test_index = X_test.index

    # run KNN for every subset with train and test splits
    for (fold_,Nsamp_),methods_ in subs.items():
        # run a subset
        if Nsamp_ not in arr_run:
            continue 
        meta_=meta[dataset_][(fold_,Nsamp_)]['metadata']
        if len(meta_.index)<Nsamp_:
            continue 
        nn_res[dataset_][(fold_,Nsamp_)]={}
        nn_res_tmp[dataset_][(fold_,Nsamp_)]={}
        for count_,(method,dist_tmp) in enumerate(methods_.items()):
            if method in ['Bray_Curtis','GUniFrac_Alpha_One', 'Robust_Aitchison']:
                nn_res[dataset_][(fold_,Nsamp_)][method]={}
                meta_=meta[dataset_][(fold_,Nsamp_)]['metadata']
                tmp_train_index=list(set(train_index)&set(meta_.index))
                tmp_test_index=list(set(test_index)&set(meta_.index))
                # encode 
                meta_['y_encode']=proc.transform(meta_[case_study[dataset_]['factor']].values).ravel()
                pcoa_tmp=pcoa(DistanceMatrix(dist_tmp)).samples[['PC1','PC2']]
                pcoa_tmp.index=meta_.index
                # Run default classifier 
                knn = KNeighborsClassifier(n_neighbors=25).fit(pcoa_tmp.loc[tmp_train_index,:], 
                                                 meta_.loc[tmp_train_index,:]['y_encode'].astype(int).ravel())
                acc_ = accuracy_score(meta_.loc[tmp_test_index,:]['y_encode'].ravel(),
                                       knn.predict(pcoa_tmp.loc[tmp_test_index,:]).astype(int))
                nn_res[dataset_][(fold_,Nsamp_)][method]['R^{2}'] = acc_
                nn_res_tmp[dataset_][(fold_,Nsamp_)]=pd.DataFrame(nn_res[dataset_][(fold_,Nsamp_)])


    both_nn[dataset_]=pd.concat(nn_res_tmp[dataset_])

make_safe_dir(subsample_results)
both_perm_res['AG'].to_csv('subsample_results/AG_types_of_plants_fstat.csv')
both_perm_res['Sponges'].to_csv('subsample_results/Sponges_health_status_fstat.csv')
both_perm_res['Sleep_Apnea'].to_csv('subsample_results/Sleep_Apnea_exposure_type_fstat.csv')

both_nn['AG'].to_csv('subsample_results/AG_types_of_plants_classifier.csv')
both_nn['Sponges'].to_csv('subsample_results/Sponges_health_status_classifier.csv')
both_nn['Sleep_Apnea'].to_csv('subsample_results/Sleep_Apnea_exposure_type_classifier.csv')

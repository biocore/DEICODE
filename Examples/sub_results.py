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
    for (fold_,Nsamp_),methods_ in subs.items():
        perm_res[dataset_][(fold_,Nsamp_)]={}
        perm_res_tmp[dataset_][(fold_,Nsamp_)]={}
        for method,dist_tmp in methods_.items():
            perm_res[dataset_][(fold_,Nsamp_)][method]={}
            dist_tmp=DistanceMatrix(dist_tmp)
            perm_tmp=permanova(dist_tmp,meta[dataset_][(fold_,Nsamp_)]['metadata'][case_study[dataset_]['factor']].values)
            perm_res[dataset_][(fold_,Nsamp_)][method]['test statistic']=perm_tmp['test statistic']
            perm_res[dataset_][(fold_,Nsamp_)][method]['p-value']=perm_tmp['p-value']
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
    for (fold_,Nsamp_),methods_ in subs.items():
        nn_res[dataset_][(fold_,Nsamp_)]={}
        nn_res_tmp[dataset_][(fold_,Nsamp_)]={}
        for method,dist_tmp in methods_.items():
            nn_res[dataset_][(fold_,Nsamp_)][method]={}
            meta_=meta[dataset_][(fold_,Nsamp_)]['metadata']
            
            # use actual loadings directly
            if method=='Robust_Aitchison':
                pcoa_tmp=ordinations[dataset_][(fold_,Nsamp_)]['RPCA_Samples']
                rename_cols={i-1:'PC'+str(i) for i in range(1,len(pcoa_tmp.columns)+1)}
                pcoa_tmp = pcoa_tmp.rename(columns=rename_cols)
                pcoa_tmp = pcoa_tmp.reindex(index=meta_.index)
            # otherwise just use the distances
            else:
                pcoa_tmp=pcoa(DistanceMatrix(dist_tmp)).samples[['PC1','PC2','PC3']]
                pcoa_tmp.index=meta_.index

            # see classieifer accuracy
            X=np.array(pcoa_tmp)
            y=meta_[case_study[dataset_]['factor']].values
            #prepare data
            X = StandardScaler().fit_transform(X)
            # prepare for classifier 
            proc=preprocessing.LabelBinarizer()
            y=proc.fit_transform(y)
            X_train, X_test, y_train, y_test = train_test_split(X, y.ravel(), test_size=.2,stratify=y, random_state=42)
            knn = KNeighborsClassifier().fit(X_train,y_train)
            acc_ = accuracy_score(y_test,knn.predict(X_test))
            #acc_train = knn.score(X_train, y_train)
            #acc_test = knn.score(X_test, y_test)            
            #nn_res[dataset_][(fold_,Nsamp_)][method]['R_sq_train'] = acc_train
            nn_res[dataset_][(fold_,Nsamp_)][method]['R_sq_test'] = acc_
            nn_res_tmp[dataset_][(fold_,Nsamp_)]=pd.DataFrame(nn_res[dataset_][(fold_,Nsamp_)])
            
    both_nn[dataset_]=pd.concat(nn_res_tmp[dataset_])

make_safe_dir(subsample_results)
both_perm_res['AG'].to_csv('subsample_results/AG_types_of_plants_fstat.csv')
both_perm_res['Sponges'].to_csv('subsample_results/Sponges_health_status_fstat.csv')
both_perm_res['Sleep_Apnea'].to_csv('subsample_results/Sleep_Apnea_exposure_type_fstat.csv')

both_nn['AG'].to_csv('subsample_results/AG_types_of_plants_classifier.csv')
both_nn['Sponges'].to_csv('subsample_results/Sponges_health_status_classifier.csv')
both_nn['Sleep_Apnea'].to_csv('subsample_results/Sleep_Apnea_exposure_type_classifier.csv')

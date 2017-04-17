#!/usr/bin/env python



from __future__ import division
#parsing command line
import argparse
import os
#pcoa
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from warnings import warn
from scipy.spatial.distance import pdist, squareform
#PCA
from sklearn.decomposition import PCA
#machine leanring
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import mpl_toolkits.mplot3d.axes3d as p3
from sklearn.multiclass import OneVsOneClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from skbio.stats.composition import clr
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
#visual
from Impute_vis import PCA_niche
import seaborn as sns
import pylab
import matplotlib
import brewer2mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#other util
import sys
from gneiss.util import match
import scipy
import pandas as pd
import numpy as np
from scipy import stats, optimize
from biom import load_table
import operator
import copy
import warnings

__author__ = 'Cameron_Martino, cameronmartino at gmail dot com'
parser = argparse.ArgumentParser()
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')

required.add_argument('-i','--Input_OTU', help='Path to .biom table i.e. home/usr/input/otu.biom (taxa should be included in metadata)',required=True)
required.add_argument('-m','--map', help='Path to Qiime style metadata i.e. home/usr/input/map.txt',required=True)
required.add_argument('-o','--output',help='Output directory', required=True)
optional.add_argument('-l',"--low_rank_method",type=str,default="SoftImpute",help='Specify a low rank method to use (default is SoftImpute) (options = NNM (Nuclear Norm Minimization), SoftImpute,IterativeSVD, MatrixFactorization)', required=False)
optional.add_argument("-d", "--decompit",type=int,default=100,help="How many iterations to complete in decomposition (default=100) (options = any integer or None)", required=False)
optional.add_argument("-b", "--bactnum",type=int,default=20,help="Number of bacteria to extract from PCA axis (default=12) (options = any integer)", required=False)
optional.add_argument("-c", "--classnum",type=int,default=3,help="Number of highest scoring classifiers to use in analysis (default=3) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
optional.add_argument("-t", "--taxause",type=str,default='species',help="What level of taxonomy to extract from PCA axis (default=genus) (options = phylum, class, order, family, genus, species or None if you do not have incorporated taxaonomy)", required=False)
optional.add_argument("-s", "--mapstart",type=int,default=0,help="What column to start analysis on in mapping file, (i.e. skipping barcode sequences) (default=0) (options = any integer greater than or equal to 0 and less than the umber of columns in the mapping file)", required=False)
optional.add_argument("-f", "--feature",type=bool,default=False,help="Set to False if you would like to turn on feature selection (default=True, warning: on large datastes this will signficantly slow down run time)", required=False)

optional.add_argument("-a", "--minsample",type=int,default=0,help="Minimum number of counts a sample needs (min=0)", required=False)
optional.add_argument("-u", "--minOTU",type=int,default=0,help="Minimum number of counts a OTU needs (min=0)", required=False)

args = parser.parse_args()

print("\n\nInput biom file is: %s"%(args.Input_OTU))
print("Input metadata is: %s"%(args.map))
print("Output directory is: %s"%(args.output))
print("iterations: %i"%(args.decompit))
print("Number of bacteria to extract from PCA axis: %i"%(args.bactnum))
print("Number of highest scoring classifiers to use: %i"%(args.classnum))
print("level of taxonomy to extract from PCA axis: %s"%(args.taxause))

in_biom=args.Input_OTU
out=args.output
map_file=args.map
lr_method=args.low_rank_method
mins=args.minsample
minotu=args.minOTU


try:
    from fancyimpute import BiScaler, NuclearNormMinimization, SoftImpute, IterativeSVD, MatrixFactorization
except:
    print("MatrixFactorization not available on this platform, please choose another method.")
    from fancyimpute import BiScaler, KNN, NuclearNormMinimization, SoftImpute, IterativeSVD
    if lr_method=="MatrixFactorization":
        sys.exit('Please choose another method, MatrixFactorization will not run')

iteration_used=args.decompit
bactnum_for_pca=args.bactnum
classnum_to_analy=args.classnum
taxause_name=args.taxause
mapstart_num=args.mapstart
select_features=args.feature

try:
    os.stat(out)
except:
    print("specified output folder does not exsist: making it now")
    os.mkdir(out)
if taxause_name == "phylum":
    txlvl=1
elif taxause_name == "class":
    txlvl=2
elif taxause_name == "order":
    txlvl=3
elif taxause_name == "family":
    txlvl=4
elif taxause_name == "genus":
    txlvl=5
elif taxause_name == "species":
    txlvl=6
elif taxause_name == "none" or taxause_name == "None":
    txlvl=42

def convert_biom_to_pandas(table): # covert biom
    feature_table = pd.DataFrame(np.array(table.matrix_data.todense()).T,index=table.ids(axis='sample'),columns=table.ids(axis='observation'))
    feature_ids = table.ids(axis='observation')
    mapping = {i: table.metadata(id=i, axis='observation')['taxonomy'] for i in feature_ids}
    for key, value in mapping.items():
        nvalue = ';'.join(value[1:])
        mapping.update({key:nvalue})
    taxonomy = pd.DataFrame(mapping, index=['taxonomy']).T
    return feature_table, taxonomy

def pw_distances(counts, ids=None, metric="braycurtis"):
    
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError("Number of rows in counts must be equal to number of provided ""ids.")
    
    distances = pdist(counts, metric)
    return DistanceMatrix(squareform(distances, force='tomatrix', checks=False), ids)

########################### metadata classification ##########################################################
print('\n Importing Metadata for analysis \n')
mappingdf= pd.read_table('%s'%map_file, index_col=0,low_memory=False)
mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)
if min(mappingdf.shape)<=0:
    sys.exit('Import error from mapping or metadata (less than two samples or features): please check your metadata is tab delimited format')
print('\n Done')

############################# import otu information ###################################################
print('\n Importing .biom table for analysis \n')
try:
    filename=in_biom.split('/')[-1]
except:
    filename=in_biom
if filename.split('.')[-1]=="biom":
    #BIOM
    #load table
    total_number_seq_sample=0
    total_number_seq_features=0
    table = load_table('%s'%in_biom)
    read_filter1 = lambda val, id_, md: sum(val) > total_number_seq_sample
    read_filter2 = lambda val, id_, md: sum(val) > total_number_seq_features
    table.filter(read_filter1, axis='sample')
    table.filter(read_filter2, axis='observation')
    otu, taxonomy = convert_biom_to_pandas(table)
    otu=otu.T
    otu=otu.replace(np.nan,0, regex=True)
    #add taxa names
    taxa_names=list(taxonomy['taxonomy'])
elif filename.split('.')[-1]=="csv" or filename.split('.')[-1]=="tsv" or filename.split('.')[-1]=="txt":
    #csv
    otu=pd.read_table('%s'%in_biom, index_col=0)
    taxa_names=list(otu.index.values)
    otu=otu.replace(np.nan,0, regex=True)
    if min(otu.shape)<=1:
        sys.exit('Import error less than two samples or features: please check that your data is tab delimited or in biom file format')
else:
    sys.exit('Import error: please check that your data is one of the following file formats (.csv,.biom,.txt,.tsv)')
#add unque taxa names for pca/machine leanring (save taxa name for later)
tax_index=[]
otus_index=[]
for q in range(0,len(otu.index.values)):
    otus_index.append("OTU_%s"%str(q))
otu['new_index']=otus_index
otu = otu.set_index('new_index')
#### match and save data #####
otu, mappingdf = match(otu.T, mappingdf)
otu=otu.T
#set taxa names ro level specified by user
if txlvl==42:
    tax_index=otus_index
else:
    for t in taxa_names:
        tax_index.append(";".join(t.split(";")[:txlvl]))
#remove otus with sum to zero after matching files
otu=otu.loc[(otu.sum(axis=1) != 0)]
# save data, names and classifiers from data frame
index = otu.index.values.tolist()
data = otu.as_matrix()
ids = otu.columns.values.tolist()
ids = list(map(str, ids))
# process taxa names
tax_index_new=[]
for cho in index:
    tax_index_new.append(tax_index[int(cho.split("_")[1])])
tax_index=tax_index_new
tax_index_new=[]
#encode pre preoccessing from mapping
samplenames = mappingdf.index.values.tolist()
samplenames = map(str, samplenames)
encoded_mapping={} #save coded and uncoded as dict
encoded_mappingdf=mappingdf.copy() # encoded dataframe
le = preprocessing.LabelEncoder() # encoder prepreocessing
classifiers_meta=mappingdf.columns.values.tolist() # classifier names
for metatmp in classifiers_meta[mapstart_num:]: # run each classifier
    le.fit(list(set(list(mappingdf[metatmp]))))
    encoded = le.transform(list(mappingdf[metatmp]))
    not_encoded = le.inverse_transform(encoded)
    encoded_mapping[metatmp]=[encoded,not_encoded]
    encoded_mappingdf[metatmp]=encoded #encoded dataframe

#size
print("Number of samples %i"%(int((data.shape[1]))))
print("Number of OTUs %i"%(int((data.shape[0]))))
print("\n Done")


############################# Low-Rank Matrix Imputations ###################################################
print("\n Running Low-Rank Matrix Imputation \n")

otum=data.copy() # make copy for imputation
otum=otum.astype(np.float64)
#test
otum[otum == 0] = np.nan #make unknown nan
if lr_method=='MatrixFactorization':
    print(" Running MatrixFactorization")
    low_rank_matrix=MatrixFactorization(min_value=0).complete(otum)
elif lr_method=='IterativeSVD':
    print(" Running Iterative SVD")
    low_rank_matrix=IterativeSVD(max_iters=iteration_used,convergence_threshold=0.000001,min_value=0).complete(otum)
elif lr_method=="NNM":
    print(" Running Nuclear Norm Minimization")
    low_rank_matrix=NuclearNormMinimization(min_value=0).complete(otum)
else:
    print(" Running Soft Impute")
    low_rank_matrix=SoftImpute(max_rank=min(otum.shape),max_iters=iteration_used,convergence_threshold=0.000001,min_value=0,max_value=(np.amax(otum))).complete(otum)
print('\nDone')

############################# SVM , based on data composition and determine best classifier ###################################################
print('\nTesting Cummultive Cumulative Explained Variance for WPCA \n')
# cumulative explained variance from weighted PCA
pca_model=PCA(n_components=3) #PCA
X_reduced_var = pca_model.fit_transform(low_rank_matrix.copy()) #transform
pccompdf = pd.DataFrame(pca_model.components_,columns=otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
var_exp=pca_model.explained_variance_ratio_
cum_var_exp = np.cumsum(pca_model.explained_variance_ratio_)
print("The Explained Variance By PCA Axis is: ")
print(var_exp)
print("The Cumulative Explained Variance By PCA is: ")
print(cum_var_exp)

############################# machine leanring ###################################################
sv={} # save scores for each classifier
split_tmp=mapstart_num
rng = np.random.RandomState(42)
for reduced_this in classifiers_meta[mapstart_num:]:

    #remove all unknown samples for metadata column and make sure it matches
    low_rank_matrix_matched,mappingdf_matched = match(pd.DataFrame(low_rank_matrix,index,ids).T,mappingdf[~mappingdf[reduced_this].isin(['Unknown'])]) #filter data for unknowns
    mappingdf_matched_encoded,mappingdf_matched = match(encoded_mappingdf,mappingdf_matched)
    low_rank_matrix_matched=low_rank_matrix_matched.as_matrix()
    
    # split data
    X_train, X_test, y_train_all, y_test_all = train_test_split(low_rank_matrix_matched.copy(),np.array(mappingdf_matched_encoded.as_matrix()),test_size=0.2,random_state=0)
    #start machine learning
    print('\nRunning Support vector machines \n')
    metatmp=reduced_this
    check_tmp=0
    for check_occurance in list(set(encoded_mapping[metatmp][1])):
        if all(isinstance(item, str) for item in encoded_mapping[metatmp][1]) and check_tmp==0:
            if list(encoded_mapping[metatmp][1]).count(check_occurance) <=1:
                check_tmp+=1
                print("    Warning: Skipping Catagory: %s contains labels that occurs only once, will cause spurious results."%(str(metatmp)))
                continue

    y_train=y_train_all.T[split_tmp]
    y_test=y_test_all.T[split_tmp]
    split_tmp+=1
    
    print("   Running Classifier: %s with %i classes"%(str(metatmp),len(set(encoded_mapping[metatmp][0]))))
    
    if len(set(encoded_mapping[metatmp][0]))<=1: # can not learn classifiers with one label
        print("    Warning: Skipping Catagory: %s,  Catagory must have more than one label!"%str(metatmp))
        continue

    if all(isinstance(item, str) for item in encoded_mapping[metatmp][1]):
        print("    Random Forest Regression")
        if select_features == True:
            clfr = RandomForestClassifier(random_state=rng)
            selector_clfr=RFE(clfr)
            selector_clfr.fit(X_train, y_train)
            sv[metatmp] = selector_clfr.score(X_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))
        else:
            clfr = RandomForestClassifier(random_state=rng)
            clfr.fit(X_train, y_train)
            sv[metatmp] = clfr.score(X_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))

    else: # qauntity
        print("    Random Forest Regression")
        if select_features == True:
            clfr = RandomForestRegressor(random_state=rng)
            selector_clfr=RFE(clfr)
            selector_clfr.fit(X_train, y_train)
            sv[metatmp] = selector_clfr.score(X_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))
        else:
            clfr = RandomForestRegressor(random_state=rng)
            clfr.fit(X_train, y_train)
            sv[metatmp] = clfr.score(X_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))

#Convert dict to dataframe and choose colors
print('\nSaving Classifier Scores and Visualizations for each classifier Based on SVM R-Squared \n')
scores = pd.DataFrame(list(sv.items()))
scores=scores.set_index(scores[0])
scores = scores.drop([0], 1)
scores.columns = ['Scores']
scores.sort_values(['Scores'], ascending = [False], inplace = True)
scores.to_csv('%s/metadata_scores.csv'%(out), sep='\t')
mybest_classer_list = scores.index.values.tolist()

#plot scores
new_names=[]
for get_s in list(scores.index.values):
    new_names.append(get_s+' (n='+str(len(mappingdf[~mappingdf[get_s].isin(['Unknown'])][get_s]))+')'+' (labels='+str(len(list(set(mappingdf[~mappingdf[get_s].isin(['Unknown'])][get_s]))))+')')
scores.index=new_names
fig, (ax1) = plt.subplots(ncols=1, nrows=1)
scores.columns=['Matrix Completion (RF)']
scores.sort_values(['Matrix Completion (RF)'], ascending = [True], inplace = True)
#rename a few
scores.plot(kind='barh',title='Mean Cross-Validation Scores (Lauber $et \, al.$)',xlim=(0,1),ax=ax1)
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1],loc=2,prop={'size':12}, bbox_to_anchor=(1.0, 1.0))
ax1.set_ylabel('')
plt.savefig('%s/machine_leanring_scores.png'%(out),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')
plt.close('all')

### clr for plot ###
low_rank_matrix_plot=clr(low_rank_matrix+1)

####### visual (clr-PCA,PCoA) #####
print('\nSaving Visualization\n')
for bestclassifier in mybest_classer_list[:classnum_to_analy]:

    #PCoA#
    
    #str
    if all(isinstance(item, str) for item in encoded_mapping[metatmp][1]):
        
        X_reduced_pcoa=pcoa(DistanceMatrix(pdist(data.T,'braycurtis'),ids)).samples[['PC1','PC2']].as_matrix()
        ids=['PC1','PC2']
        spPCAplot_pcoa=pd.DataFrame(X_reduced_pcoa,low_rank_matrix_plot_matched.index.values,ids)
        spPCAplot_pcoa['labels']=encoded_mapping[bestclassifier][1].tolist()
        sns.swarmplot(x="PC1", y="PC2", data=spPCAplot_pcoa, hue="labels", size=10,ax=ax1)
        ax1.set_xlabel('$PC-1$')
        ax1.set_ylabel('$PC-2$')
        ax1.legend(loc=2,prop={'size':22},bbox_to_anchor=(1, 1))
        plt.savefig('%s/%s_pca.png'%(out,bestclassifier),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')
        plt.close('all')

    else: #cont
        
        Y=encoded_mapping[bestclassifier][1].tolist()
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        bc_dm=pw_distances(data.T.copy(), ids)
        ord_results=pcoa(bc_dm)
        X_reduced = ord_results.samples.as_matrix()
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCoA")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("$PC-1$", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("$PC-2$", fontsize=10)
        ax.set_yticklabels([])
        plt.colorbar(p)
        fig.set_tight_layout(True)
        plt.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)
        plt.close('all')


    #PCA#

    #str
    if all(isinstance(item, str) for item in encoded_mapping[metatmp][1]):
        
        pca_model=PCA(n_components=2)
        X_reduced2 = pca_model.fit_transform(low_rank_matrix_plot.T.as_matrix())
        ids=['PC1','PC2']
        spPCAplot_pca=pd.DataFrame(X_reduced2,low_rank_matrix_plot_matched.index.values,ids)
        spPCAplot_pca['labels']=encoded_mapping[bestclassifier][1].tolist()
        sns.swarmplot(x="PC1", y="PC2", data=spPCAplot_pca, hue="labels", size=10,ax=ax1)
        ax1.set_xlabel('$PC-1$')
        ax1.set_ylabel('$PC-2$')
        ax1.legend(loc=2,prop={'size':22},bbox_to_anchor=(1, 1))
        plt.savefig('%s/%s_pca.png'%(out,bestclassifier),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')
        plt.close('all')

    else: #cont
        
        Y=encoded_mapping[bestclassifier][1].tolist()
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        reduce_clf=PCA(n_components=2)
        X_reduced = reduce_clf.fit_transform(low_rank_matrix_plot.T.copy())
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCA")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("$PC-1$", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("$PC-2$", fontsize=10)
        ax.set_yticklabels([])
        plt.colorbar(p)
        fig.set_tight_layout(True)
        plt.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)
        plt.close('all')

############################# Plot OTUS affecting each classifier for best axis of PCA ########################################################

    # Niche bacterial variance visualization

    PC_list=['PC-1','PC-2']
    for pc in PC_list:
        # otu orignal data
        if all(isinstance(item, str) for item in encoded_mapping[metatmp][1]):
            cont=False # false if str
        else:
            cont=True
        # Extract information from imputed PCA axis
        out_niche_linkeddf,observed_table_sfi,index_mean,index_std,highest_var_bact,pccompdf = PCA_niche.niche_visual(otu,low_rank_matrix,tax_index,bactnum_for_pca,bestclassifier,pc,mappingdf)
        out_niche_linkeddf.to_csv('%s/%s_bact_plot_%s.csv'%(out,bestclassifier,pc), sep='\t')
        pccompdf.to_csv('%s/%s_most_variance_bact_%s.csv'%(out,bestclassifier,pc), sep='\t')
        # Visualize the data
        plt = PCA_niche.plot_niche(out_niche_linkeddf,observed_table_sfi,mappingdf,encoded_mapping,bestclassifier,pc,index_mean,index_std,le,cont)
        plt.savefig('%s/%s_bacteria_extarcted_axis_%s.png'%(out,bestclassifier,pc),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')

    plt.close("all")

print('\n Done, thank you for using DEICODE \n\n')

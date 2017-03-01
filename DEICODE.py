#!/usr/bin/env python

from __future__ import division
#parsing command line
import argparse
import os
#low rank methods (new)
from wpca import WPCA, EMPCA
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
optional.add_argument('-l',"--low_rank_method",type=str,default="SoftImpute",help='Specify a low rank method to use (default is SoftImpute) (options = NNM (Nuclear Norm Minimization), SoftImpute,IterativeSVD, MatrixFactorization, can also use WPCA, or EMPCA for imputation)', required=False)
optional.add_argument("-d", "--decompit",type=int,default=200,help="How many iterations to complete in decomposition (default=100) (options = any integer or None)", required=False)
optional.add_argument("-b", "--bactnum",type=int,default=10,help="Number of bacteria to extract from PCA axis (default=12) (options = any integer)", required=False)
optional.add_argument("-c", "--classnum",type=int,default=3,help="Number of highest scoring classifiers to use in analysis (default=2) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
optional.add_argument("-t", "--taxause",type=str,default='class',help="What level of taxonomy to extract from PCA axis (default=genus) (options = phylum, class, order, family, genus, species or None if you do not have incorporated taxaonomy)", required=False)
optional.add_argument("-s", "--mapstart",type=int,default=0,help="What column to start analysis on in mapping file, (i.e. skipping barcode sequences) (default=0) (options = any integer greater than or equal to 0 and less than the umber of columns in the mapping file)", required=False)
optional.add_argument("-f", "--feature",type=bool,default=True,help="Set to False if you would like to turn off feature selection (default=True, warning: on large datastes this will signficantly slow down run time)", required=False)
optional.add_argument("-w_zero", "--w_zero",type=int,default=0,help="Zero value used for weighted PCA (default=0)", required=False)
optional.add_argument("-w_low", "--w_low",type=int,default=.00001,help="Low value used for weighted PCA (default .001)", required=False)
optional.add_argument("-w_high", "--w_high",type=int,default=1000,help="High value used for weighted PCA (default 10)", required=False)
optional.add_argument("-n", "--ncomp",type=int,default=3,help="Number of Principle Components to Use in Feature Selection (default=3)", required=False)
optional.add_argument("-fm", "--fmethod",type=str,default="WPCA",help="Method to use for feature selection (default=WPCA)", required=False)
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
zerow=args.w_zero
minw=args.w_low
maxw=args.w_high
component=args.ncomp
feature_method=args.fmethod

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
#Mapping import
mappingdf= pd.read_table('%s'%map_file, index_col=0,low_memory=False)
mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)
#mappingdf = mappingdf[mappingdf.life_stage != 'Unknown'] #TODO add metadata filtering feature 
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
    total_number_seq_sample=0 #TODO add to input
    total_number_seq_features=0 #TODO add to input
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

#set taxa names ro level specified by user
if txlvl==42:
    tax_index=otus_index
else:
    for t in taxa_names:
        tax_index.append(";".join(t.split(";")[:txlvl]))

#### match and save data #####

otu, mappingdf = match(otu.T, mappingdf)
otu=otu.T

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
print("Number of samples %i"%(int(min(data.shape))))
print("Number of OTUs %i"%(int(max(data.shape))))
print("\n Done")


############################# Low-Rank Matrix Imputations ###################################################

print("\n Running Low-Rank Matrix Imputation \n")
if lr_method=="WPCA" or lr_method=="EMPCA": # Impute by WPCA or EMPCA
    
    # higher frequency is determined as having less noise
    weight=data.T.copy()
    weight = (weight - weight.min(axis=0)) / (weight.max(axis=0) - weight.min(axis=0)) * (maxw - zerow) + zerow

    if lr_method=="EMPCA":
        print(" Running EMPCA")
        imputem=EMPCA(n_components=3)
        low_rank_matrix = imputem.fit_reconstruct(data.T.copy(),weight).T
    else:
        print(" Running WPCA")
        imputem=WPCA(n_components=3)
        low_rank_matrix = imputem.fit_reconstruct(data.T.copy(),weight).T

else:
    
    otum=data.copy() # make copy for imputation
    otum=otum.astype(np.float64)
    #test
    otum[otum == 0] = np.nan #make unknown nan
    if lr_method=='MatrixFactorization':
        print(" Running MatrixFactorization")
        low_rank_matrix=MatrixFactorization().complete(otum)
    elif lr_method=='IterativeSVD':
        print(" Running Iterative SVD")
        low_rank_matrix=IterativeSVD(rank=min(otum.shape),max_iters=iteration_used,convergence_threshold=0.000001).complete(otum)
    elif lr_method=="NNM":
        print(" Running Nuclear Norm Minimization")
        low_rank_matrix=NuclearNormMinimization().complete(otum)
    else:
        print(" Running Soft Impute")
        low_rank_matrix=SoftImpute(max_rank=min(otum.shape),max_iters=iteration_used,convergence_threshold=0.000001,min_value=0,max_value=(np.amax(otum))).complete(otum)

print('\nDone')

############################# SVM , based on data composition and determine best classifier ###################################################


print('\nTesting Cummultive Cumulative Explained Variance for WPCA \n')

if lr_method=="EMPCA" or lr_method=="WPCA":

    X_reduced_var = imputem.fit_transform(low_rank_matrix.copy(),weight) #transform
    pccompdf = pd.DataFrame(imputem.components_,columns=otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
    var_exp=imputem.explained_variance_ratio_
    cum_var_exp = np.cumsum(imputem.explained_variance_ratio_)
    print("The Explained Variance By PC Axis is: ")
    print(var_exp)
    print("The Cumulative Explained Variance is: ")
    print(cum_var_exp)

else:
    # cumulative explained variance from weighted PCA
    weight=low_rank_matrix.copy()
    weight = (weight - weight.min(axis=0)) / (weight.max(axis=0) - weight.min(axis=0)) * (maxw - zerow) + zerow
    pca_model=WPCA(n_components=3) #PCA
    X_reduced_var = pca_model.fit_transform(low_rank_matrix.copy(),weight) #transform
    pccompdf = pd.DataFrame(pca_model.components_,columns=otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
    var_exp=pca_model.explained_variance_ratio_
    cum_var_exp = np.cumsum(pca_model.explained_variance_ratio_)
    print("The Explained Variance By PC Axis is: ")
    print(var_exp)
    print("The Cumulative Explained Variance is: ")
    print(cum_var_exp)


# feature selection method
if select_features == True:
    if feature_method=="WPCA":
        feature_clf = WPCA(n_components=component)
    else:
        feature_clf = EMPCA(n_components=component)

# split data
X_train, X_test, y_train_all, y_test_all = train_test_split(low_rank_matrix.copy().T,np.array(encoded_mappingdf.as_matrix()),test_size=0.2,random_state=0)

#feature selection
if select_features == True:
    print("\n Running Weighted PCA for Feature Selection and Dimensionality Reduction \n")
    weights=X_train.copy()
    weights = (weights - weights.min(axis=0)) / (weights.max(axis=0) - weights.min(axis=0)) * (maxw - zerow) + zerow
    feature_clf.fit(X_train,weights)
    X_t_train = feature_clf.transform(X_train)
    X_t_test = feature_clf.transform(X_test)

#start machine learning
print('\nRunning Support vector machines \n')

sv={} # save scores for each classifier
split_tmp=mapstart_num
for metatmp in classifiers_meta[mapstart_num:]: # run each classifier
    
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
    
    elif len(set(encoded_mapping[metatmp][0]))==2: # if only two possible classifications use support vector classfier
        
        print("    Running Support Vector Classifier")

        if select_features == True:
            clfb = svm.SVC(random_state=0)
            clfb.fit(X_t_train, y_train)
            sv[metatmp] = clfb.score(X_t_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))
        else:
            clfb = svm.SVC(random_state=0)
            clfb.fit(X_train, y_train)
            sv[metatmp] = clfb.score(X_test, y_test)

    elif len(set(encoded_mapping[metatmp][0]))>2 and all(isinstance(item, str) for item in encoded_mapping[metatmp][1]) and len(set(encoded_mapping[metatmp][0]))<200: # if not quantity and class is not boolian
    
        print("    Running One vs. One Linear Support Vector Classifier")

        if select_features == True:
            
            # feature selection
            clfm = OneVsOneClassifier(svm.LinearSVC(random_state=0,class_weight="balanced")) # one vs one for imblanced data sets
            clfm.fit(X_t_train, y_train)
            sv[metatmp] = clfm.score(X_t_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))
        
        else:
            clfm = OneVsOneClassifier(svm.LinearSVC(random_state=0,class_weight="balanced")) # one vs one for imblanced data sets
            clfm.fit(X_train, y_train)
            sv[metatmp] = clfm.score(X_test, y_test)

    else: # if qauntity
        
        print("    Running Support Vector Regression")
        if select_features == True:
            
            X_t_train = feature_clf.transform(X_train)
            X_t_test = feature_clf.transform(X_test)
            clfr = GridSearchCV(svm.SVR(kernel='rbf'),cv=5,param_grid={"C": [1e0, 1e1, 1e2, 1e3],"gamma": np.logspace(-2, 2, 5)}, n_jobs=-1) # grid search for optimized perams
            clfr.fit(X_t_train, y_train)
            sv[metatmp] = clfr.score(X_t_test, y_test)
            print("    Coef: %f"%(sv[metatmp]))
        else:
            #learn
            clfr = GridSearchCV(svm.SVR(kernel='rbf',gamma=0.1),cv=5,param_grid={"C": [1e0, 1e1, 1e2, 1e3],"gamma": np.logspace(-2, 2, 5)}, n_jobs=-1) # grid search for optimized perams
            clfr.fit(X_train, y_train)
            sv[metatmp] = clfr.score(X_test, y_test)


#Convert dict to dataframe and choose colors
print('\nSaving Classifier Scores and Visualizations for each classifier Based on SVM R-Squared \n')
scores = pd.DataFrame(list(sv.items()))
scores=scores.set_index(scores[0])
scores = scores.drop([0], 1)
scores.columns = ['Scores']
scores.sort_values(['Scores'], ascending = [False], inplace = True)
scores.to_csv('%s/metadata_scores.csv'%(out), sep='\t')
mybest_classer_list = scores.index.values.tolist()
bmap = brewer2mpl.get_map('Set3','qualitative',12,reverse=True)
colors = bmap.mpl_colors

print('\nSaving Visualization\n')

### visualize ###

for bestclassifier in mybest_classer_list[:classnum_to_analy]:
    
    # bray-curtis PCOA comparison to WPCA
    
    #  continuous data with color bar
    
    if len(set(encoded_mapping[bestclassifier][0])) > 20 and ( all(isinstance(item, int) for item in encoded_mapping[bestclassifier][1]) or all(isinstance(item, float) for item in encoded_mapping[bestclassifier][1])):

        Y=encoded_mapping[bestclassifier][1].tolist()
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        bc_dm=pw_distances(data.T.copy(), ids)
        ord_results=pcoa(bc_dm)
        X_reduced = ord_results.samples.as_matrix()
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCA on Origonal Matrix: colored by pH")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_yticklabels([])
        plt.colorbar(p)
        fig.set_tight_layout(True)
        plt.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)
        plt.close('all')

        # non-continuous data PCOA with legend
        
    else:
        
        Y=encoded_mapping[bestclassifier][1].tolist()
        bc_dm=pw_distances(data.T.copy(), ids)
        ord_results=pcoa(bc_dm)
        X_reduced = ord_results.samples.as_matrix()
        bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
        colors = bmap7.mpl_colors
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        scatter_proxy=[]
        k=set(Y)
        v=colors
        #qlists
        q1=[0]
        q2=[0,8]
        q3=[0,4,8]
        q4=[0,3,5,8]
        q5=[0,2,4,6,8]
        q6=[0,1,3,5,6,8]
        q7=[0,2,3,5,6,7,8]
        q8=[0,1,2,3,4,5,7,8]
        q9=[0,1,2,3,4,5,6,7,8]
        choodict={1: q1,2: q2,3: q3,4: q4,5: q5,6: q6,7: q7,8: q8,9: q9}
        #choose qlist
        for key,value in choodict.items():
            if int(key)==int((len(set(Y)))):
                for q in value:
                    scatter_proxy.append(matplotlib.lines.Line2D([0],[0], linestyle="none", c=colors[q], marker = 'o'))

        if int(len(k))>=3:
            numrows=4
        else:
            numrows=int(len(k))
        ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows, labelspacing=.1, borderaxespad=0.,prop={'size':25})
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.Set1_r,s=150)
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_yticklabels([])
        fig.set_tight_layout(True)
        fig.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)



	# plot PCA in 3D

	#discr
    if len(set(encoded_mapping[bestclassifier][1])) > 12:
        
        Y=encoded_mapping[bestclassifier][1].tolist()
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        X_reduced = feature_clf.fit_transform(low_rank_matrix.T.copy())
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCA on Origonal Matrix: colored by pH")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_yticklabels([])
        plt.colorbar(p)
        fig.set_tight_layout(True)
        plt.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)
        plt.close('all')

	#cont
    else:
        Y=encoded_mapping[bestclassifier][1].tolist()
        X_reduced = feature_clf.fit_transform(low_rank_matrix.T.copy())
        bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
        colors = bmap7.mpl_colors
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 15),sharey=False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows, labelspacing=.01, borderaxespad=0.,prop={'size':25})
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.Set1_r,s=150)
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.set_xticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_yticklabels([])
        fig.set_tight_layout(True)
        fig.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)


############################# Plot OTUS affecting each classifier for best axis of PCA ########################################################

	# Niche bacterial variance visualization

    PC_list=['PC-1','PC-2']
    for pc in PC_list:
        # otu orignal data
        if len(set(encoded_mapping[bestclassifier][1])) > 12:
            cont=True # true, false if not continous data
        else:
            cont=False
        # Extract information from imputed PCA axis
        out_niche_linkeddf,observed_table_sfi,index_mean,index_std,highest_var_bact,pccompdf = PCA_niche.niche_visual(otu,low_rank_matrix,tax_index,bactnum_for_pca,bestclassifier,pc,mappingdf,weight)
        out_niche_linkeddf.to_csv('%s/skin_plot_%s.csv'%(out,pc), sep='\t')
        pccompdf.to_csv('%s/most_variance_bact_%s.csv'%(out,pc), sep='\t')
        # Visualize the data
        plt = PCA_niche.plot_niche(out_niche_linkeddf,observed_table_sfi,mappingdf,encoded_mapping,bestclassifier,pc,index_mean,index_std,le,cont,weight)
        plt.savefig('%s/%s_bacteria_extarcted_axis_%s.png'%(out,bestclassifier,pc),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')

plt.close("all")

print('\n Done, thank you for using DEICODE \n\n')

#!/usr/bin/env python


from __future__ import division
#parsing command line
import argparse
#low rank methods (new)
from wpca import WPCA, EMPCA
from fancyimpute import BiScaler, KNN, NuclearNormMinimization, SoftImpute, IterativeSVD, MICE
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
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import silhouette_score
from sklearn.multiclass import OneVsOneClassifier
#visual
from Impute_vis import PCA_niche
import seaborn as sns
import pylab
import matplotlib
import brewer2mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#other util
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

parser = argparse.ArgumentParser(description='Multivariate analysis of 16S data through DEICODE\n\n')
parser.add_argument('-i','--Input_OTU', help='Path to .biom table i.e. home/usr/input/otu.biom (taxa should be included in metadata)',required=True)
parser.add_argument('-m','--map', help='Path to Qiime style metadata i.e. home/usr/input/map.txt',required=True)
parser.add_argument('-o','--output',help='Output directory', required=True)
parser.add_argument('-l','--low_rank_method',type=str,default="IterativeSVD",help='Specify a low rank method to use (default is IterativeSVD) (options = NNM (Nuclear Norm Minimization), SoftImpute,IterativeSVD, MICE, KNN, WPCA, or EMPCA)', required=False)
parser.add_argument("-d", "--decompit",type=int,default=1000,help="How many iterations to complete in decomposition (deafault=1000) (options = any integer or None)", required=False)
parser.add_argument("-b", "--bactnum",type=int,default=10,help="Number of bacteria to extract from PCA axis (default=12) (options = any integer)", required=False)
parser.add_argument("-c", "--classnum",type=int,default=3,help="Number of highest scoring classifiers to use in analysis (default=2) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
parser.add_argument("-t", "--taxause",type=str,default='genus',help="What level of taxonomy to extract from PCA axis (deafult=genus) (options = phylum, class, order, family, genus, species or None if you do not have incorporated taxaonomy)", required=False)
parser.add_argument("-s", "--mapstart",type=int,default=3,help="What column to start analysis on in mapping file, (i.e. skipping barcode sequences) (deafult=3) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
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
iteration_used=args.decompit
bactnum_for_pca=args.bactnum
classnum_to_analy=args.classnum
taxause_name=args.taxause
mapstart_num=args.mapstart

if taxause_name == "phylum":
    txlvl=2
elif taxause_name == "class":
    txlvl=3
elif taxause_name == "order":
    txlvl=4
elif taxause_name == "family":
    txlvl=5
elif taxause_name == "genus":
    txlvl=6
elif taxause_name == "species":
    txlvl=7
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
mappingdf= pd.read_table('%s'%map_file, index_col=0)
mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)
samplenames = mappingdf.index.values.tolist()
samplenames = map(str, samplenames)

#encode pre preoccessing from mapping
encoded_mapping={} #save coded and uncoded
le = preprocessing.LabelEncoder() # encoder prepreocessing

classifiers_meta=mappingdf.columns.values.tolist() # classifier names
for metatmp in classifiers_meta[mapstart_num:]: # run each classifier
    le.fit(list(set(list(mappingdf[metatmp]))))
    encoded = le.transform(list(mappingdf[metatmp]))
    not_encoded = le.inverse_transform(encoded)
    encoded_mapping[metatmp]=[encoded,not_encoded]

print('\n Done')

############################# import otu information ###################################################


print('\nImporting .biom table for analysis \n')


#BIOM

#load table
table = load_table('%s'%in_biom)
read_filter = lambda val, id_, md: sum(val) > 0
table.filter(read_filter, axis='sample')
table.filter(read_filter, axis='observation')
otu, taxonomy = convert_biom_to_pandas(table)
otu=otu.T
otu=otu.replace(np.nan,0, regex=True)


#add taxa names
taxa_names=list(taxonomy['taxonomy'])
tax_index=[]

#add unque taxa names for now
otus_index=[]
for q in range(len(otu.index.values)):
    otus_index.append("OTU_%s"%str(q))
otu['new_index']=otus_index
otu = otu.set_index('new_index')


#set taxa names ro level specified by user
if txlvl==42:
    # Generate "fake" OTU names for table
    tax_index=otus_index
else:
    for t in taxa_names: # taxa level split and join
        tax_index.append(";".join(t.split(";")[:txlvl]))

#check matching
otu, mappingdf = match(otu.T, mappingdf)
otu=otu.T

# save data and names from data frame
index = otu.index.values.tolist()
data = otu.as_matrix()
ids = otu.columns.values.tolist()
ids = list(map(str, ids))

print('\n Done')


############################# Low-Rank Matrix Imputations ###################################################


print('\n Running Low-Rank Matrix Imputation \n')




if lr_method=="WPCA" or lr_method=="EMPCA": # WPCA or EMPCA specified
    
    otum=data.T.copy() # make copy for imputation
    #WPCA and EMPCA
    weight = otum.copy()
    for i in range(len(otum)):
        for j in range(len(otum[i])):
            if otum[i][j]==0:
                weight[i][j]=1
            else:
                weight[i][j]=10000

    if lr_method=="EMPCA":
        print("Running EMPCA")
        low_rank_matrix = EMPCA(n_components=3).fit_reconstruct(otum,weight).T
    else:
        print("Running WPCA")
        low_rank_matrix = WPCA(n_components=3).fit_reconstruct(otum,weight).T

else:
    
    otum=data.copy() # make copy for imputation
    # Fancy Impute
    # make zero unknown
    otum=otum.astype(np.float64)
    otum[otum == 0] = np.nan #make unknown nan
    # fancy impute methods (SoftImpute,IterativeSVD, MICE, KNN)
    if lr_method=='KNN':
        print("Running KNN")
        low_rank_matrix=KNN(k=12,orientation="rows",use_argpartition=True,print_interval=100,min_value=0,max_value=(np.amax(otum)/10),normalizer=None,verbose=True).complete(otum)
    elif lr_method=='IterativeSVD':
        print("Running Iterative SVD")
        low_rank_matrix=IterativeSVD(convergence_threshold=0.00001,max_iters=iteration_used,gradual_rank_increase=False,svd_algorithm="arpack",init_fill_method="zero",min_value=1e-10,max_value=(np.amax(otum)),verbose=True).complete(otum)
    elif lr_method=='MICE':
        print("Running MICE")
        low_rank_matrix=MICE(visit_sequence='monotone',n_imputations=100,n_burn_in=40,n_pmm_neighbors=40,impute_type='col',n_nearest_columns=np.infty,init_fill_method="mean",min_value=0,max_value=(np.amax(otum)),verbose=True).complete(otum)
    elif lr_method=="NNM":
        print("Running Nuclear Norm Minimization")
        low_rank_matrix=NuclearNormMinimization(min_value=0,max_value=(np.amax(otum))).complete(otum)
    else:
        print("Running Soft Impute")
        low_rank_matrix=SoftImpute(shrinkage_value=None,convergence_threshold=0.0001,max_iters=iteration_used,n_power_iterations=1,init_fill_method="zero",min_value=0,max_value=(np.amax(otum)),normalizer=None,verbose=True).complete(otum)

print('\nDone')

############################# SVM , based on data composition and determine best classifier ###################################################

print('\n Testing Cummultive Cumulative Explained Variance for PCA \n')


X=low_rank_matrix.copy() # low rank matrix to run SVM on

pca_model=PCA(n_components=3) #PCA
X_reduced_var = pca_model.fit_transform(X) #transform
pccompdf = pd.DataFrame(pca_model.components_,columns=otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
var_exp=pca_model.explained_variance_ratio_
cum_var_exp = np.cumsum(pca_model.explained_variance_ratio_)
#plot cummulitave variance
with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(6,4))
    plt.bar(range(3), var_exp, alpha=0.5, align='center',
    label='Individual explained variance')
    plt.step(range(3), cum_var_exp, where='mid',label='cumulative explained Variance')
    plt.title(('Imputed Explained variance\n'))
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.legend(loc='best')
    plt.savefig('%s/explained_variance.png'%(out),bbox_to_anchor=(2.2, 1.0), dpi=300, bbox_inches='tight')
    plt.tight_layout()
#plot distribution
g = sns.jointplot(x="PC-1", y="PC-2", data=pccompdf, kind="kde", color="m")
g.plot_joint(plt.scatter, c="w", s=30, linewidth=1, marker="+")
g.ax_joint.collections[0].set_alpha(0)
g.set_axis_labels("$PC-1$", "$PC-2$")
g.savefig('%s/PCA_dist.png'%(out),dpi=300)

print('\nDone: PCA Evaluation Graphs in Output')

print('\nRunning Support vector machines \n')

sv={} # save scores for each classifier

for metatmp in classifiers_meta[mapstart_num:]: # run each classifier
    
    print("   Running Classifier: %s"%str(metatmp))

if len(set(encoded_mapping[metatmp][0]))<=1: # can not learn classifiers with one label
        print("   Wanring: Skipping Classifier: %s, All classifiers must have more than one label!"%str(metatmp))
        continue

    if len(set(encoded_mapping[metatmp][0]))==2: # if only two possible classifications use support vector classfier

        Y=encoded_mapping[metatmp][0]
        X_train, X_test, y_train, y_test = train_test_split(X.T, Y, test_size=0.2, random_state=0)
        pca = PCA(n_components=3)
        pca.fit(X_train)
        X_t_train = pca.transform(X_train)
        X_t_test = pca.transform(X_test)
        clf = svm.SVC(class_weight='balanced',random_state=0)
        clf.fit(X_t_train, y_train)
        sv[metatmp] = clf.score(X_t_test, y_test)

    if len(set(encoded_mapping[metatmp][0]))>2 and 10>len(set(encoded_mapping[metatmp][0])): # if not continuous but more than two use linear support vector classifier

        Y=encoded_mapping[metatmp][0]
        X_train, X_test, y_train, y_test = train_test_split(X.T, Y, test_size=0.2, random_state=0)
        pca = PCA(n_components=3)
        pca.fit(X_train)
        X_t_train = pca.transform(X_train)
        X_t_test = pca.transform(X_test)
        clf = OneVsOneClassifier(svm.LinearSVC(random_state=0,class_weight='balanced'))
        clf.fit(X_t_train, y_train)
        sv[metatmp] = clf.score(X_t_test, y_test)

    else: # if continuous classifier use regression
	
        Y=encoded_mapping[metatmp][0]
        X_train, X_test, y_train, y_test = train_test_split(X.T, Y, test_size=0.2, random_state=0)
        pca = PCA(n_components=3)
        pca.fit(X_train)
        X_t_train = pca.transform(X_train)
        X_t_test = pca.transform(X_test)
        clf = svm.SVR()
        clf.fit(X_t_train, y_train)
        sv[metatmp] = clf.score(X_t_test, y_test)


############ Convert dict to dataframe and choose colors ##################################


print('\n Saving Classifier Scores and Visualizations for each classifier Based on SVM R-Squared \n')

scores = pd.DataFrame(list(sv.items()))
scores=scores.set_index(scores[0])
scores = scores.drop([0], 1)
scores.columns = ['Scores']
bmap = brewer2mpl.get_map('Set3','qualitative',12,reverse=True)
colors = bmap.mpl_colors
scores.sort_values(['Scores'], ascending = [False], inplace = True)
scores.to_csv('%s/metadata_scores.csv'%out, sep='\t')
mybest_classer_list = scores.index.values.tolist()

rint('\nDone: Ranked Classifier Scores in Output')


for bestclassifier in mybest_classer_list[:classnum_to_analy]:
    
    # bray-curtis PCOA comparison to PCA
    
    #  continuous data with color bar
    if len(set(encoded_mapping[bestclassifier][0])) > 12:
        Y=encoded_mapping[bestclassifier][1].tolist()
        X=data.T
        fig = plt.figure(2, figsize=(10, 8))
        ax = Axes3D(fig,elev=-150, azim=110)
        bc_dm=pw_distances(X, ids)
        ord_results=pcoa(bc_dm)
        X_reduced = ord_results.samples.as_matrix()
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCA on Origonal Matrix: colored by pH")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("PC3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        plt.colorbar(p)
        plt.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)
        plt.close('all')

        # non-continuous data PCOA with legend
        
    else:
        Y=encoded_mapping[bestclassifier][1].tolist()
        X =data.T
        bc_dm=pw_distances(X, ids)
        ord_results=pcoa(bc_dm)
        X_reduced = ord_results.samples.as_matrix()
        bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
        colors = bmap7.mpl_colors
        fig = plt.figure()
        ax = Axes3D(fig, elev=-135, azim=230)
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

        ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows,fontsize = 'large', labelspacing=.1, borderaxespad=0.)
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.Set1_r,s=150)
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("PC3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        fig.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)



	# plot PCA in 3D

	#discr
    if len(set(encoded_mapping[bestclassifier][1])) > 12:
        Y=encoded_mapping[bestclassifier][1].tolist()
        X =low_rank_matrix.T
        fig = plt.figure(2, figsize=(10, 8))
        ax = Axes3D(fig,elev=-150, azim=110)
        X_reduced = PCA(n_components=3).fit_transform(X)
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
        ax.set_title("PCA on Origonal Matrix: colored by pH")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("PC3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        plt.colorbar(p)
        plt.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)
        plt.close('all')

	#cont
    else:
        Y=encoded_mapping[bestclassifier][1].tolist()
        X =low_rank_matrix.T
        X_reduced = PCA(n_components=3).fit_transform(X)
        bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
        colors = bmap7.mpl_colors
        fig = plt.figure()
        ax = Axes3D(fig, elev=-135, azim=230)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows,fontsize = 'large', labelspacing=.01, borderaxespad=0.)
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.Set1_r,s=150)
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("PC1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("PC2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("Pc3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        fig.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)


############################# Plot OTUS affecting each classifier for best axis of PCA ########################################################

	# Niche bacterial variance visualization

    PC_list=['PC-1','PC-2']
    for pc in PC_list:
        # otu orignal data
        imputed_in=low_rank_matrix.copy() # imputed data to use
        if len(set(encoded_mapping[bestclassifier][1])) > 12:
            cont=True # true, false if not continous data
        else:
            cont=False
        # Extract information from imputed PCA axis
        out_niche_linkeddf,observed_table_sfi,index_mean,index_std,highest_var_bact,pccompdf = PCA_niche.niche_visual(otu,low_rank_matrix,tax_index,bactnum_for_pca,bestclassifier,pc,mappingdf)
        out_niche_linkeddf.to_csv('%s/skin_plot_%s.csv'%(out,pc), sep='\t')
        pccompdf.to_csv('%s/most_variance_bact_%s.csv'%(out,pc), sep='\t')
        # Visualize the data
        plt = PCA_niche.plot_niche(out_niche_linkeddf,observed_table_sfi,mappingdf,encoded_mapping,bestclassifier,pc,index_mean,index_std,le,cont)
        plt.savefig('%s/%s_bacteria_extarcted_axis_%s.png'%(out,bestclassifier,pc),bbox_to_anchor=(1.0, 1.0), dpi=300, bbox_inches='tight')

plt.close("all")

print('\n Done, thank you for using DEICODE \n\n')

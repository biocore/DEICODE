#!/usr/bin/env python

from __future__ import division
import argparse
import scipy
import pandas as pd
import numpy as np
import itertools
from random import shuffle
import random
from scipy import stats, optimize
from skbio import DistanceMatrix
from skbio.stats.ordination import PCoA
from skbio.stats.distance import pwmantel
from skbio.diversity.beta import pw_distances
from skbio.stats.distance import bioenv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skbio.stats.distance import ANOSIM
import brewer2mpl
from skbio.stats.spatial import procrustes 
import seaborn as sns
import pylab
import matplotlib
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn import svm
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import mpl_toolkits.mplot3d.axes3d as p3
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import silhouette_score
from sklearn.ensemble import RandomForestClassifier
from biom import load_table
import operator

__author__ = 'Cameron_Martino'

parser = argparse.ArgumentParser(description='Multivariate analysis of 16S data through DEICODE\n\n')
parser.add_argument('-i','--input_dir', help='path to .biom table i.e. home/usr/input/otu.biom',required=True)
parser.add_argument('-m','--map', help='path to Qiime style metadata i.e. home/usr/input/map.txt',required=True)
parser.add_argument('-o','--output',help='Output directory', required=True)
parser.add_argument("-d", "--decompit",type=int,default=48,help="How many iterations to complete in decomposition (deafault=48) (options = any integer)", required=False)
parser.add_argument("-b", "--bactnum",type=int,default=12,help="Number of bacteria to extract from PCA axis (default=12) (options = any integer)", required=False)
parser.add_argument("-c", "--classnum",type=int,default=6,help="Number of highest scoring classifiers to use in analysis (default=2) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
parser.add_argument("-t", "--taxause",type=str,default='taxonomy',help="What level of taxonomy to extract from PCA axis (deafult=family) (options = phylum, class, order, family, genus, species)", required=False)
parser.add_argument("-s", "--mapstart",type=int,default=3,help="What column to start analysis on in mapping file, (i.e. skipping barcode sequences) (deafult=3) (options = any integer greater than 1 and less than the umber of columns in the mapping file)", required=False)
args = parser.parse_args()
print "\n\nInput biom file is: %s"%(args.input_dir)
print "Input metadata is: %s"%(args.map)
print "Output directory is: %s"%(args.output)
print "iterations: %i"%(args.decompit)
print "Number of bacteria to extract from PCA axis: %i"%(args.bactnum)
print "Number of highest scoring classifiers to use: %i"%(args.classnum)
print "level of taxonomy to extract from PCA axis: %s"%(args.taxause)
in_biom=args.input_dir
out=args.output
map_file=args.map
iteration_used=args.decompit
bactnum_for_pca=args.bactnum
classnum_to_analy=args.classnum
taxause_name=args.taxause
mapstart_num=args.mapstart


class R_pca:
    
	""" 

	COD written by Papanicolaou Alex, RPCA, (2011), GitHub repository, https://github.com/apapanico/RPCA 

	"""

	def __init__(self, D, mu=None, lmbda=None):
		self.D = D
		self.S = np.zeros(self.D.shape)
		self.Y = np.zeros(self.D.shape)


		if mu:
			self.mu = mu
		else:
			self.mu = np.prod(self.D.shape) / (4 * self.norm_p(self.D, 2))

		self.mu_inv = 1 / self.mu

		if lmbda:
			self.lmbda = lmbda
		else:
			self.lmbda = 1 / np.sqrt(np.max(self.D.shape))

	@staticmethod
	def norm_p(M, p):
		return np.sum(np.power(M, p))

	@staticmethod
	def shrink(M, tau):
		return np.sign(M) * np.maximum((np.abs(M) - tau), np.zeros(M.shape))

	def svd_threshold(self, M, tau):
		U, S, V = np.linalg.svd(M, full_matrices=False)
		return np.dot(U, np.dot(np.diag(self.shrink(S, tau)), V))

	def fit(self, tol=None, max_iter=1000, iter_print=100):
		iter = 0
		err = np.Inf
		Sk = self.S
		Yk = self.Y
		Lk = np.zeros(self.D.shape)

		if tol:
			_tol = tol
		else:
			_tol = 1E-7 * self.norm_p(np.abs(self.D), 2)

		while (err > _tol) and iter < max_iter:
			Lk = self.svd_threshold(self.D - Sk + self.mu_inv * Yk, self.mu_inv)
			Sk = self.shrink(self.D - Lk + (self.mu_inv * Yk), self.mu_inv * self.lmbda)
			Yk = Yk + self.mu * (self.D - Lk - Sk)
			err = self.norm_p(np.abs(self.D - Lk - Sk), 2)
			iter += 1
			if (iter % iter_print) == 0 or iter == 1 or iter > max_iter or err <= _tol:
				print 'iteration: {0}, error: {1}'.format(iter, err)

		self.L = Lk
		self.S = Sk
		return Lk, Sk

def convert_biom_to_pandas(table):

    """ 
    Biom to pandas dataframe code:

    Jamie Morton, gneiss, (2016), GitHub repository, https://github.com/biocore/gneiss


    Unpacks biom table into two pandas dataframes.

    The first dataframe will contain the count information for 
    features and samples. The second datafram will contain taxonomy 
    information for all of the OTUs.

    Parameters
    ----------
    table : biom.Table

    Returns
    -------
    pd.DataFrame
    Contingency table of counts where samples correspond 
    to rows and columns correspond to features (i.e. OTUs)
    pd.DataFrame
    A mapping of OTU names to taxonomic ids
    """

    feature_table = pd.DataFrame(np.array(table.matrix_data.todense()).T,index=table.ids(axis='sample'),columns=table.ids(axis='observation'))
    feature_ids = table.ids(axis='observation')
    mapping = {i: table.metadata(id=i, axis='observation')['taxonomy'] for i in feature_ids}
    for key, value in mapping.iteritems():
        nvalue = ';'.join(value[1:])
        mapping.update({key:nvalue})
    # modify below as necessary.  
    # There are typically 7 levels of taxonomy.
    taxonomy = pd.DataFrame(mapping, index=['taxonomy']).T
    return feature_table, taxonomy




########################### metadata classification ##########################################################

print '\nImporting map.txt for analysis \n'
print '\nDone'
df= pd.read_table('%s'%map_file, index_col=0)
label_save=[]
#save names
map_names= pd.read_table('%s'%map_file, index_col=0)


# set sting based netadata to numerical classifiers

classifier_names = df.columns.values.tolist()
imdex_list = df.index.values.tolist()
for nm in classifier_names:

    if df[nm].dtype == np.dtype(object) or df[nm].dtype == np.dtype(bool) :
        
        # get set of strings in column

        dfList = df[nm].tolist()
        dflistset = set(dfList)
        dfsetdict={}
        q=0
        
        # assign each string in set a integer classifier

        for tmpclassify in dflistset:
            dfsetdict[tmpclassify]=q
            q+=1
        label_save.append(dfsetdict)
        # re name data according assigned values

       
        for in_tmp in imdex_list:
            tmp_rename=df[nm][str(in_tmp)]
            for key, value in dfsetdict.iteritems():
                if tmp_rename == key:
                    df.set_value(str(in_tmp), nm, value)
                    
     
samplenames = df.index.values.tolist()
samplenames = map(str, samplenames)

############################# import otu information ###################################################


print '\nImporting .biom table for analysis \n'
print '\nDone'

table = load_table('%s'%in_biom)
otu, taxonomy = convert_biom_to_pandas(table)
otu=otu.T

#add taxa names

taxonomy['fg'] = taxonomy[taxause_name]
result = pd.concat([otu, taxonomy], axis=1, join='inner')
otu2 = result.set_index('fg')
#otu2 = otu2.drop(['kingdom','phylum','class','order','family','genus','species'], 1)
otu=otu2

# reoder according to classifier Y

otu=otu.T
otu=otu.reindex(samplenames)
otu=otu.T
otu_check=otu   
df = df.fillna(0)
otu =otu.fillna(0)

# save data and names from data frame

index = otu.index.values.tolist()
data = otu.as_matrix()
ids = otu.columns.values.tolist()
ids = map(str, ids)

#check mapping again 

df=df.reindex(ids)

############################# Decompose ###################################################

print '\nRunning Decomposition \n'

rpcadata=np.array(data)
rpca = R_pca(rpcadata)
L, S = rpca.fit(max_iter=iteration_used, iter_print=int(iteration_used/4))
low_rank_matrix = L
sparse_matrix = S
save_lowrankdf = pd.DataFrame(low_rank_matrix,index=index, columns=ids)
save_sparsedf = pd.DataFrame(sparse_matrix, columns=ids)
print '\nDone'
############################# SVM , based on data composition and determine best classifier ###################################################

# fill

print '\nRunning Support vector machines \n' 
print '\nDone'
classifiers_meta=df.columns.values.tolist()
X =low_rank_matrix.T
sv={}

for metatmp in classifiers_meta[mapstart_num:]: 

	if len(set(metatmp))==2:

		Y=df[metatmp].tolist()
		Y = map(int, Y)
		X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, Y, test_size=0.2, random_state=0)
		pca = PCA(n_components=3)
		pca.fit(X_train)
		X_t_train = pca.transform(X_train)
		X_t_test = pca.transform(X_test)
		clf = svm.SVC()
		clf.fit(X_t_train, y_train)
		sv[metatmp] = clf.score(X_t_test, y_test)

	if len(set(metatmp))>2 and 4>len(set(metatmp)):


		Y=df[metatmp].tolist()
		Y = map(int, Y)
		X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, Y, test_size=0.2, random_state=0)
		pca = PCA(n_components=3)
		pca.fit(X_train)
		X_t_train = pca.transform(X_train)
		X_t_test = pca.transform(X_test)
		clf = svm.LinearSVC(random_state=29)
		clf.fit(X_t_train, y_train)
		sv[metatmp] = clf.score(X_t_test, y_test)

	else:
	
		Y=df[metatmp].tolist()
		Y = map(int, Y)
		X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, Y, test_size=0.2, random_state=0)
		pca = PCA(n_components=3)
		pca.fit(X_train)
		X_t_train = pca.transform(X_train)
		X_t_test = pca.transform(X_test)
		clf = svm.SVR()
		clf.fit(X_t_train, y_train)
		sv[metatmp] = clf.score(X_t_test, y_test)


############ Convert dict to dataframe and choose colors ##################################


print '\nSaving otuput files and Visualizations for each classifier \n' 

scores = pd.DataFrame(sv.items())
scores=scores.set_index(scores[0])
scores = scores.drop([0], 1)
scores.columns = ['Scores']
bmap = brewer2mpl.get_map('Set3','qualitative',12,reverse=True)
colors = bmap.mpl_colors
scores.sort_values(['Scores'], ascending = [False], inplace = True)
scores.to_csv('%s/metadata_scores.csv'%out, sep='\t')
mybest_classer_list = scores.index.values.tolist()

#save best classifiers 

for bestclassifier in mybest_classer_list[:classnum_to_analy]:



	############################# Print OTUS affecting each classifier for best axis of PCA ########################################################
    
    #PCOA
    if len(set(df[bestclassifier].tolist())) > 12:
        Y=df[bestclassifier].tolist()
        X =rpcadata.T
        fig = plt.figure(2, figsize=(10, 8))
        ax = Axes3D(fig,elev=-150, azim=110)
        bc_dm=pw_distances(X, ids, "braycurtis")
        ord_results=PCoA(bc_dm).scores()
        X_reduced = ord_results.site.T
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.cool,s=200)
        p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.cool,s=200)
        ax.set_title("PCA on Origonal Matrix: colored by pH")
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("F1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("F2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("F3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        plt.colorbar(p)
        plt.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)
        plt.close('all')
    

    
        #cont
    else:
        Y=df[bestclassifier].tolist()
        Y2=map_names[bestclassifier].tolist()
        X =rpcadata.T
        bc_dm=pw_distances(X, ids, "braycurtis")
        ord_results=PCoA(bc_dm).scores()
        X_reduced = ord_results.site.T
        bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
        colors = bmap7.mpl_colors
        fig = plt.figure()
        ax = Axes3D(fig, elev=-120, azim=230)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        scatter_proxy=[]
        k, v = set(Y2), colors
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
        for key,value in choodict.iteritems():
            if key==int((len(set(Y2)))):
                for q in value:
                    scatter_proxy.append(matplotlib.lines.Line2D([0],[0], linestyle="none", c=v[q], marker = 'o'))

        if int(len(k))>=3:
            numrows=4
        else:
            numrows=int(len(k))

        ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows,fontsize = 'large', labelspacing=.1, borderaxespad=0.)
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.Set1_r,s=150)
        ax.set_axis_bgcolor('white')
        ax.set_xlabel("F1", fontsize=10)
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel("F2", fontsize=10)
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel("F3", fontsize=10)
        ax.w_zaxis.set_ticklabels([])
        fig.savefig('%s/%s_pcoa.png'%(out,bestclassifier), dpi=300)



	# plot PCA in 3D

	#discr

	if len(set(df[bestclassifier].tolist())) > 12:
		Y=df[bestclassifier].tolist()
		X =low_rank_matrix.T
		fig = plt.figure(2, figsize=(10, 8))
		ax = Axes3D(fig,elev=-150, azim=110)
		X_reduced = PCA(n_components=3).fit_transform(X)
		ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.cool,s=200)
		p=ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.cool,s=200)
		ax.set_title("PCA on Origonal Matrix: colored by pH")
		ax.set_axis_bgcolor('white')
		ax.set_xlabel("F1", fontsize=10)
		ax.w_xaxis.set_ticklabels([])
		ax.set_ylabel("F2", fontsize=10)
		ax.w_yaxis.set_ticklabels([])
		ax.set_zlabel("F3", fontsize=10)
		ax.w_zaxis.set_ticklabels([])
		plt.colorbar(p)
		plt.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)
		plt.close('all')




	#cont 
	else:
		Y=df[bestclassifier].tolist()
		Y2=map_names[bestclassifier].tolist()
		X =low_rank_matrix.T
		X_reduced = PCA(n_components=3).fit_transform(X)
		bmap7 = brewer2mpl.get_map('Set1','qualitative',int(len(set(Y2))),reverse=True)
		colors = bmap7.mpl_colors
		fig = plt.figure()
		ax = Axes3D(fig, elev=-120, azim=230)
		handles, labels = ax.get_legend_handles_labels()
		ax.legend(handles, labels)
		ax.legend(scatter_proxy, k, numpoints = 1,bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numrows,fontsize = 'large', labelspacing=.01, borderaxespad=0.)
		ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,cmap=plt.cm.Set1_r,s=150)
		ax.set_axis_bgcolor('white')
		ax.set_xlabel("F1", fontsize=10)
		ax.w_xaxis.set_ticklabels([])
		ax.set_ylabel("F2", fontsize=10)
		ax.w_yaxis.set_ticklabels([])
		ax.set_zlabel("F3", fontsize=10)
		ax.w_zaxis.set_ticklabels([])
		fig.savefig('%s/%s_pca.png'%(out,bestclassifier), dpi=300)
		

	# extract OTUs from each axis 

	# PCA

	pca = PCA(n_components=3)
	pca.fit_transform(save_lowrankdf.T)

	# Dump components relations with features:

	df2 = pd.DataFrame(pca.components_,columns=save_lowrankdf.T.columns,index = ['PC-1','PC-2','PC-3']).T
	df2.abs()
	PC_list=['PC-1','PC-2','PC-3']
	for pc in PC_list:


		df2.sort_values([pc], ascending = [False], inplace = True)
		bact_list = df2.index.values.tolist()
		check='start'

		bactandotu=otu


		for bact in bact_list[:bactnum_for_pca]:
		    if check=='start':
			def3=bactandotu.T[bact]
			def3.sum(axis=1)
			def2=def3.T[-1:]
			def2=def2.T
			def2['%s_labels'%bestclassifier] = df[bestclassifier]
			def2=def2.groupby('%s_labels'%bestclassifier, as_index=False).mean()
			def2 = def2.set_index('%s_labels'%bestclassifier)
			check='not-start'
		    else:
			def3tmp=bactandotu.T[bact]
			def3tmp.sum(axis=1)
			def2tmp=def3tmp.T[-1:]
			def2tmp=def2tmp.T
			def2tmp['%s_labels'%bestclassifier] = df[bestclassifier]
			def2tmp=def2tmp.groupby('%s_labels'%bestclassifier, as_index=False).mean()
			def2tmp = def2tmp.set_index('%s_labels'%bestclassifier)
			def2 = pd.concat([def2, def2tmp], axis=1, join='inner')

		for ln in def2.columns.values.tolist():
			if ln=='f__':
				def2=def2.drop(ln, 1)
		    
		def2=def2.T
		def2=def2.groupby(def2.index).sum()
		def2=def2.T 


		def2.to_csv('%s/%s_variance_axis_%s.csv'%(out,bestclassifier,pc), sep='\t')


		fig, (ax1) = plt.subplots(ncols=1, nrows=1, figsize=(36, 10))
		bmap = brewer2mpl.get_map('Set3','qualitative',12,reverse=True)
		colors = bmap.mpl_colors
		def2.plot(kind='bar', stacked=False,color=colors,legend=True, width=.80, figsize=(13, 10),fontsize=26,ax=ax1)
		if int(int(len(set(map_names[bestclassifier].tolist())))/4) == 0:
			numcol=1
		else:
			numcol=int(int(len(set(map_names[bestclassifier].tolist())))/4)
		ax1.legend(bbox_to_anchor=(0., .89, 1., .102), loc=3, ncol=numcol,fontsize = 'large', labelspacing=.2, borderaxespad=0.)
		ax1.set_xticklabels(ax1.xaxis.get_majorticklabels(), rotation=360)
		ax1.tick_params(axis='y', labelsize=26)
		ax1.set_axis_bgcolor('white')
		ax1.set_ylabel("OTU Count Averaged by Sub Group",{'size':'26'})
		ax1.set_xlabel("%s"%bestclassifier,{'size':'26'})
		plt.savefig('%s/%s_bacteria_extarcted_axis_%s.png'%(out,bestclassifier,pc), dpi=300)

plt.close("all")
print '\nDone, thank you for using DEICODE \n\n'

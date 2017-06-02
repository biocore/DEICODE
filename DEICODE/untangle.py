#utils
from __future__ import division
import pandas as pd
import numpy as np
from scipy import stats, optimize
import sys
import itertools
#ML
from sklearn import preprocessing
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.cluster.bicluster import SpectralCoclustering
#PCOA
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
#transforms 
from skbio.stats.composition import clr,centralize
#PCA 
from sklearn.decomposition import PCA
#ploting
import seaborn as sns
import matplotlib.pyplot as plt
#completion
from fancyimpute import SoftImpute
#DEICODE utils
from . import fetch

def _scatter(ax1,reduced,catvis,legendadd=True):

    '''
    input: OTU table and underlying rank structure 
    
    output: sorted dataframe with orginal index and column labels, dataframe with assigned niche labels 
    
    '''
    
    for ((key, grp)) in reduced.groupby(catvis):
        ax1.scatter(grp['PC1'], grp['PC2'], color=next(ax1._get_lines.prop_cycler)['color'], label=key, s=50)
    ax1.set_xlabel('$PC-1$')
    ax1.set_ylabel('$PC-2$')
    ax1.set_xticks([])
    ax1.set_yticks([])
    if legendadd==True:
        ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))

    return 

def biplot(data,r):
    
    '''
    input: OTU table and underlying rank structure 
    
    output: two dataframes both coclustered, first one with orginal row and column labels, second with assigned row and column labels 
    
    '''
    
    #generate model 
    model = SpectralCoclustering(n_clusters=r, random_state=0)
    model.fit(data.copy())
    #sort 
    data['_sort_']=model.row_labels_
    data.sort_values(by='_sort_',inplace=True)
    data.drop('_sort_', axis=1, inplace=True)
    data=data.T
    data['_sort_']=model.column_labels_
    data.sort_values(by='_sort_',inplace=True)
    data.drop('_sort_', axis=1, inplace=True)
    
    return data.T,pd.DataFrame(data.as_matrix().T,columns=list(np.sort(model.column_labels_)),index=list(np.sort(model.row_labels_)))


def compositional(df):
    
    '''
    input: dataframe with raw counts (sprase)
    
    output: relative abundance table to (i.e. each column sums to 100% )
    
    '''

    comp_count=df.copy()
    for n in comp_count.columns.values.tolist():
        if n != 'bact_id':
            wtemp = (100/(comp_count["%s"%n].sum()))
            comp_count['%s'%n] = comp_count['%s'%n] *wtemp
    return comp_count

def encode_mapping(mapdf):
    
    '''
    
    input: mapping dataframe
    
    output: encoded mapping dataframe
    
    '''
    
    mapdf[list(mapdf.columns.values)] = mapdf[list(mapdf.columns.values)].astype(str)
    return mapdf.apply(preprocessing.LabelEncoder().fit_transform)

def complete_matrix(data,iteration,minval=0):
    
    '''
    input: otu table data (sprase)
    
    output: completed otu table (dense)
    
    '''
    
    otum=data.copy() # make copy for imputation
    otum=otum.astype(np.float64)
    otum[otum == 0] = np.nan #make unknown nan
    return SoftImpute(max_rank=min(otum.shape),max_iters=iteration,convergence_threshold=0.00001,min_value=minval,max_value=(np.amax(otum)),verbose=False).complete(otum)

def scatter_plot(reduced,mapping,catvis):
    
    '''
        input: data dataframe, mapping file dataframe, mapping file column to plot
        
        output: scatter plot for mapping data
        
        '''
    
    if all(isinstance(item, str) for item in list(mapping.T[catvis])) or all(isinstance(item, bool) for item in list(mapping.T[catvis])):
        fig, (ax1) = plt.subplots(ncols=1, nrows=1, figsize=(8, 6),sharey=False)
        for ((key, grp)) in reduced.groupby(catvis):
            plt.scatter(grp['PC1'], grp['PC2'], color=next(ax1._get_lines.prop_cycler)['color'], label=key, s=50)
        ax1.set_xlabel('$PC-1$')
        ax1.set_ylabel('$PC-2$')
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.legend(loc=2,prop={'size':22},bbox_to_anchor=(1.0, 1.0))
        return fig
    else:
        fig, (ax1) = plt.subplots(ncols=1, nrows=1, figsize=(8, 6),sharey=False)
        reduced.plot.scatter(x='PC1', y='PC2', c=str(catvis), s=50,ax=ax1);
        ax1.set_xlabel('$PC-1$')
        ax1.set_ylabel('$PC-2$')
        ax1.set_xticks([])
        ax1.set_yticks([])
        return fig

def feature_vis(otulearn,mappingdf,importance,catv,taxa):
    
    groupbothotum=compositional(otulearn).T[importance]
    groupbothotum,taxa=fetch.matchtable(groupbothotum,taxa.T)
    groupbothotum.columns=taxa.T['taxonomy']
    groupbothotum=groupbothotum.T
    
    new_taxa=[]
    #lowest taxaclassification + (phylum)
    for nm in groupbothotum.index.values.tolist():
        #find lowest taxa level
        q=0
        while q<(len(nm.split(';'))-1):
            lw=nm.split(';')[q]
            q+=1
            if len(nm.split(';')[q])==3 or '1'in nm.split(';')[q] or '2' in nm.split(';')[q]:
                break
        if '[' in lw.split('__')[1]:
            lw=lw.split('__')[1]
            lw=lw[1:-1]
        else:
            lw=lw.split('__')[1]
        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
    groupbothotum.index=new_taxa
    groupbothotum=groupbothotum.T

    if all(isinstance(item, str) for item in list(mappingdf.T[catv])):
        groupbothotum['group']=mappingdf.T[[catv]]
        dfs=[]
        name=[]
        for region, df_region in groupbothotum.groupby('group'):
            dftmp=df_region
            dftmp.drop('group', axis=1, inplace=True)
            dfs.append(dftmp)
            name.append(region)
        fig, axs = plt.subplots(ncols=1, nrows=len(dfs),sharey=False)
        for count, df in enumerate(dfs):
            axt=axs[count]
            label=name[count]
            df.plot(kind='bar',stacked=True, alpha = 0.92,rot=90,colormap="Set1",legend=False,sharey=False,fontsize=15,ax=axt)
            axt.set_ylabel('$Relative Abundance$ (%)',fontsize=12)
            axt.set_xlabel('Samples Within Group %s'%(label),fontsize=12)
            axt.set_xticks([])
    else:
        groupbothotum.index=list(mappingdf.T[catv])
        groupbothotum.index=list(map(int,groupbothotum.index))
        fig, ax1 = plt.subplots(ncols=1, nrows=1,sharey=True)
        groupbothotum.groupby(by=groupbothotum.index).mean().plot(kind='area', alpha = 0.92,lw=6,rot=45,colormap="Set1",legend=True,sharey=True,fontsize=15,ax=ax1)
        ax1.set_ylabel('$Relative Abundance$ (%)',fontsize=18)
        ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))

    if all(isinstance(item, str) for item in list(mappingdf.T[catv])):
        axt.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 2.0))

    return fig
        
def features_ml(otulearn,mapdf,catv,complete=True,iteration=100):
    
    rng = np.random.RandomState(42)
    
    try:
        mapdf[[catv]]
    except:
        sys.exit('Value given is not a mapping catagory')
    
    #mapping 
    mapdfch=mapdf.copy()
    mapdfen=encode_mapping(mapdf)

    
    if all(isinstance(item, str) for item in mapdfch[catv]) or all(isinstance(item, bool) for item in mapdfch[catv]):
        clfr = RandomForestClassifier(random_state=rng)
        clfr.fit(otulearn.copy().T.as_matrix(), mapdf[catv].values.ravel())
        importance = clfr.feature_importances_
        importance = pd.DataFrame(importance, index=otulearn.index,columns=["Importance"])
        importance["Std"] = np.std([tree.feature_importances_ for tree in clfr.estimators_], axis=0)
        importance=importance.sort_values(by=['Importance'],ascending=False)
        return importance
         
    else: 
        clfr = RandomForestRegressor(random_state=rng)
        clfr.fit(otulearn.copy().T.as_matrix(), mapdf[catv].values.ravel())
        importance = clfr.feature_importances_
        importance = pd.DataFrame(importance, index=otulearn.index,columns=["Importance"])
        importance["Std"] = np.std([tree.feature_importances_ for tree in clfr.estimators_], axis=0)
        importance=importance.sort_values(by=['Importance'],ascending=False)

    return importance
        

def get_confusion(otulearn,mapdf,catv,splits=0.2):

    '''
    input: Otu table and mapping file in dataframe and catagory in mapping
    
    output: normalized confusion matrix (dataframe)
    
    '''
    #match
    otulearn,mapdf=fetch.matchtable(otulearn,mapdf.T)
    rng = np.random.RandomState(42)
    X_train, X_test, y_train, y_test = train_test_split(otulearn.copy().T.as_matrix(),mapdf.T[catv].values.ravel(),test_size=splits,random_state=0)      
    # Compute confusion matrix
    clfr = RandomForestClassifier(random_state=rng)
    clfr.fit(X_train, y_train)
    y_pred=clfr.predict(X_test)
    cnf_matrix = confusion_matrix(y_test, y_pred)
    cnf_matrix = cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis]

    return pd.DataFrame(cnf_matrix,index=list(set(mapdf.T[catv].values)),columns=list(set(mapdf.T[catv].values))).astype(float)


    
def machine_learning(otulearn,mapdf,complete=True,single=False,single_cat=[''],iteration=100,mean_count=10,addtofilter=[]):
    
    '''
    input: matched otu tables and mapping data 
    
    output: Random Forest machine learning per catagory in mapping data
    
    '''
    
    #generate filter list for bad labels from mapping (i.e. unkown)
    if not addtofilter:
        filterlist=['Unknown']
    else:
        filterlist=['Unknown']+addtofilter
    
    #complete matrix 
    if complete==True:
        otulearn=pd.DataFrame(complete_matrix(otulearn.as_matrix(),iteration),columns=otulearn.columns,index=otulearn.index)

    
    #intitals 
    sv={} # save scores for each classifier
    split_tmp=0
    rng = np.random.RandomState(0)
    
    if single==True:
        split_tmp=[list(mapdf.columns).index(s) for s in single_cat]
        X_train, X_test, y_train_all, y_test_all = train_test_split(otulearn.as_matrix().copy().T,np.array(mapdf.as_matrix()),test_size=0.2,random_state=0)
        y_train=y_train_all.T[split_tmp]
        y_test=y_test_all.T[split_tmp]
        single_cat
        clfr = RandomForestClassifier(random_state=rng)
        clfr.fit(X_train, y_train)
        sv[reduced_this] = clfr.score(X_test, y_test)
        #Convert dict to dataframe and choose colors
        scores = pd.DataFrame(list(sv.items()))
        scores=scores.set_index(scores[0])
        scores = scores.drop([0], 1)
        scores.columns = ['Scores']
        scores.sort_values(['Scores'], ascending = [False], inplace = True)
        #add info to names
        new_names=[]
        for get_s in list(scores.index.values):
            new_names.append(get_s+' (n='+str(len(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))+')'+' (labels='+str(len(list(set(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))))+')')
        scores.index=new_names
        return scores,otulearn
                
    #run ml
    for reduced_this in list(mapdf.columns):
        
        #mapping
        tmp_otutabledf,tmpmapdfch=fetch.matchtable(otulearn,mapdf[~mapdf[reduced_this].isin(filterlist)].T)
        tmpmapdfch=tmpmapdfch.T
        mapdfch=tmpmapdfch.copy()
        tmpmapdf_encode=encode_mapping(tmpmapdfch)
        
        q=0 #take mean of n iterations
        crsvl=[] #store them
        while q!=mean_count:
            
            # split
            X_train, X_test, y_train_all, y_test_all = train_test_split(tmp_otutabledf.as_matrix().copy().T,np.array(tmpmapdf_encode.as_matrix()),test_size=0.2,random_state=0)

            # check quality of classifier
            if len(set(mapdf[reduced_this]))<=1: # can not learn classifiers with one label
                break
            
            #set reduced labels 
            y_train=y_train_all.T[split_tmp]
            y_test=y_test_all.T[split_tmp]

            # learn 
            if all(isinstance(item, str) for item in mapdf[reduced_this]) or all(isinstance(item, bool) for item in mapdf[reduced_this]):
                clfr = RandomForestClassifier(random_state=rng)
                clfr.fit(X_train, y_train)
                crsvl.append(clfr.score(X_test, y_test))
            else: 
                clfr = RandomForestRegressor(random_state=rng)
                clfr.fit(X_train, y_train)
                crsvl.append(clfr.score(X_test, y_test)) 
            
            q+=1
        if not crsvl:
            continue
        else:
            sv[reduced_this]=np.mean(crsvl)
        split_tmp+=1
    
    #Convert dict to dataframe and choose colors
    scores = pd.DataFrame(list(sv.items()))
    scores=scores.set_index(scores[0])
    scores = scores.drop([0], 1)
    scores.columns = ['Scores']
    scores.sort_values(['Scores'], ascending = [False], inplace = True)
    #add info to names
    new_names=[]
    for get_s in list(scores.index.values):
        new_names.append(get_s+' (n='+str(len(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))+')'+' (labels='+str(len(list(set(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))))+')')
    scores.index=new_names
    
    return scores,otulearn

def reduce_plot(otudf,mapping,catvis,method='clrPCA',iters=200,weights_extract=True,reduce_only=False,min_val=1e-15):
    
    '''
        input: matched otu tables and mapping data and a column of the mapping file for dimentionality reduction
        
        output: 2D scatter plot for specified column (All) and variant features for each axis as dataframes (not for PCoA)
        
        '''
    
    if method=='PCoA':
        pcaplot=pcoa(DistanceMatrix(pdist(otudf.as_matrix().T,'braycurtis'),list(otudf.columns))).samples[['PC1','PC2','PC3']]
        pcaplot[catvis]=list(mapping.T[catvis])
        if reduce_only==True:
            return pcaplot
        else:
            return scatter_plot(pcaplot,mapping,catvis)
    
    elif method=='PCA':
        pca_model=PCA(n_components=3)
        reduced=pca_model.fit_transform(otudf.as_matrix().T)
        pcaplot=pd.DataFrame(reduced,columns=['PC1','PC2','PC3'],index=otudf.columns)
        pcaplot[catvis]=list(mapping.T[catvis])
        OTUweights=pd.DataFrame(pca_model.components_.T,columns=['PC1','PC2','PC3'],index=otudf.index)
        if reduce_only==True:
            return pcaplot
        if weights_extract==True:
            return scatter_plot(pcaplot,mapping,catvis),OTUweights.sort(['PC1'],ascending=False),OTUweights.sort(['PC2'],ascending=False)
        else:
            return scatter_plot(pcaplot,mapping,catvis)
    else:
        pca_model=PCA(n_components=3)
        reduced=pca_model.fit_transform(clr(complete_matrix(otudf.as_matrix(),iteration=iters,minval=min_val)).T)
        pcaplot=pd.DataFrame(reduced,columns=['PC1','PC2','PC3'],index=otudf.columns)
        pcaplot[catvis]=list(mapping.T[catvis])
        OTUweights=pd.DataFrame(pca_model.components_.T,columns=['PC1','PC2','PC3'],index=otudf.index)
        if reduce_only==True:
            return pcaplot
        if weights_extract==True:
            return scatter_plot(pcaplot,mapping,catvis),OTUweights.sort(['PC1'],ascending=False),OTUweights.sort(['PC2'],ascending=False)
        else:
            return scatter_plot(pcaplot,mapping,catvis)

#utils
from __future__ import division
import pandas as pd
import numpy as np
import sys
#ML
from sklearn import preprocessing
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split # TODO remove after single fix
from sklearn.cluster.bicluster import SpectralCoclustering
from sklearn.model_selection import ShuffleSplit
#ploting
import matplotlib.pyplot as plt
#completion/
from fancyimpute import SoftImpute
# utils
from gneiss.util import match



def complete_matrix(data,iteration,minval):
    
    
    """
        
    Replace all zeros with small non-zero values. Using a collaborative filtering based matrix completion method.
    A replacement for adding only small values to data before performing transforms.
    Also useful for removing sparsity constraints when performing downstream analysis.
    
    ----------
    
    data: array_like a matrix of counts
    rows = Features (i.e. OTUs, metabolites)
    columns = Samples
    
    iteration: float, optional
    The number of convex iterations to optomize the solution
    If iteration is not specified, then the default iteration is 100. Which redcues to a satisfactory error threshold.
    
    
    minval: float, optional
    A small number to be used to replace zeros
    If minval is not specified, then the default minval is 1e-3. Worked well in practice with compositional transforms.
    
    
    
    Returns
    -------
    numpy.ndarray, np.float64
    A completely dense matrix of counts
    
    
    Raises
    ------
    ValueError
    Raises an error if input is a pandas dataframe and not a numpy array
    `ValueError: Lengths must match to compare`.
    
    
    Notes
    -----
    Assumes a low-rank underlying matrix, this means it performs poorly in gradient like tables. Future high-rank completion methods can overcome this.
    
    
    References
    ----------
    .. [1] Rubinsteyn A, Feldman S. 2016. fancyimpute: Version 0.0.16.
    .. [2] Mazumder R, Hastie T, Tibshirani R. 2010. Spectral Regularization Algorithms for Learning Large Incomplete Matrices. J Mach Learn Res 11:2287–2322.
    .. [3] Pending Publication; Martino and Morton
    
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.impute import complete
    >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
    >>> complete(X)
    array([[ 0.2       ,  0.4       ,  0.4       ,  0.001     ],
    [ 0.23683603,  0.5       ,  0.5       ,  0.00118418]])
    
    """
    
    otum=data.copy().astype(np.float64) # make copy for imputation, check type
    otum[otum == 0] = np.nan # make all previously zero values unknown (i.e. Nan)
    # return imputed matrix through Fancy Impute's Soft Impute function
    return SoftImpute(max_rank=min(otum.shape),max_iters=iteration,convergence_threshold=0.00001,min_value=minval,max_value=(np.amax(otum)),verbose=False).complete(otum)

def biplot(data,r,time=False):
    
    '''
    input: OTU table in dataframe, underlying rank structure r , (if time == True columns are not sorted, returns in input time)
    
    output: Two dataframes both coclustered, first one with orginal row and column labels, second with assigned row and column labels
    
    '''
    #generate model 
    model = SpectralCoclustering(n_clusters=r, random_state=0)
    model.fit(data.copy())
    #sort 
    data['_sort_']=model.row_labels_
    data.sort_values(by='_sort_',inplace=True)
    data.drop('_sort_', axis=1, inplace=True)
    data=data.T
    if time==False:
        data['_sort_']=model.column_labels_
        data.sort_values(by='_sort_',inplace=True)
        data.drop('_sort_', axis=1, inplace=True)
        return data.T,pd.DataFrame(data.as_matrix().T,columns=list(np.sort(model.column_labels_)),index=list(np.sort(model.row_labels_)))
    elif time==True:
        return data.T,pd.DataFrame(data.as_matrix().T,columns=list(model.column_labels_),index=list(np.sort(model.row_labels_)))

def relative_abund_(df):
    
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


def feature_vis(otulearn,mappingdf,importance,catv,taxa):
    
    groupbothotum=relative_abund_(otulearn.T).T[importance].T
    groupbothotum,taxa=match(groupbothotum,taxa)
    groupbothotum.index=taxa['taxonomy']
    new_taxa=[]
    # lowest taxaclassification + (phylum)
    for nm in groupbothotum.index.values.tolist():
        # find lowest taxa level
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
        new_taxa.append(lw+' ('+nm.split(';')[1].split('__')[1]+')')
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
        fig, axs = plt.subplots(ncols=1, nrows=len(dfs),sharey=True)
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
    
    
    try:
        mapdf[[catv]]
    except:
        sys.exit('Value given is not a mapping catagory')
    
    #mapping 
    mapdfch=mapdf.copy()
    mapdfen=encode_mapping(mapdf)

    if all(isinstance(item, str) for item in mapdfch[catv]) or all(isinstance(item, bool) for item in mapdfch[catv]):
        clfr = RandomForestClassifier(random_state=0)
        clfr.fit(otulearn.copy().T.as_matrix(), mapdf[catv].values.ravel())
        importance = clfr.feature_importances_
        importance = pd.DataFrame(importance, index=otulearn.index,columns=["Importance"])
        importance["Std"] = np.std([tree.feature_importances_ for tree in clfr.estimators_], axis=0)
        importance=importance.sort_values(by=['Importance'],ascending=False)
        return importance
         
    else: 
        clfr = RandomForestRegressor(random_state=0)
        clfr.fit(otulearn.copy().T.as_matrix(), mapdf[catv].values.ravel())
        importance = clfr.feature_importances_
        importance = pd.DataFrame(importance, index=otulearn.index,columns=["Importance"])
        importance["Std"] = np.std([tree.feature_importances_ for tree in clfr.estimators_], axis=0)
        importance=importance.sort_values(by=['Importance'],ascending=False)

    return importance

def machine_learning(otulearn,mapdf,complete=True,single=False,single_cat=[''],nmin=1e-3,test_split=0.2,iteration=100,mean_count=10,addtofilter=[]):
    
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
        otulearn=pd.DataFrame(complete_matrix(otulearn.as_matrix(),iteration,minval=nmin),columns=otulearn.columns,index=otulearn.index)
    
    #intitals 
    sv={} # save scores for each classifier
    split_tmp=0
    
    if single==True:
    
        #filter all columns for filterlist
        for current_ in single_cat:
            mapdf=mapdf[~mapdf[current_].isin(filterlist)]
        mapdf=mapdf[single_cat]
        #mapping
        tmp_otutabledf,tmpmapdfch=match(otulearn.T,mapdf)
        tmp_otutabledf=tmp_otutabledf.T
        mapdfch=tmpmapdfch.copy()
        tmpmapdf_encode=encode_mapping(tmpmapdfch)
        
        crsvl=[] #store them
        sss = ShuffleSplit(n_splits=mean_count, test_size=test_split, random_state=0)
        sss.get_n_splits(tmp_otutabledf.as_matrix().copy().T, np.array(tmpmapdf_encode.as_matrix()))
        for train_index, test_index in sss.split(tmp_otutabledf.as_matrix().copy().T, np.array(tmpmapdf_encode.as_matrix())):
            X_train, X_test = tmp_otutabledf.as_matrix().copy().T[train_index], tmp_otutabledf.as_matrix().copy().T[test_index]
            y_train, y_test = np.array(tmpmapdf_encode.as_matrix())[train_index].T, np.array(tmpmapdf_encode.as_matrix())[test_index].T
            
            # learn
            if all(isinstance(item, str) for item in mapdf[single_cat]) or all(isinstance(item, bool) for item in mapdf[single_cat]):
                clfr = RandomForestClassifier(random_state=0, bootstrap=False)
                clfr.fit(X_train, y_train.T)
                crsvl.append(clfr.score(X_test, y_test.T))
            else:
                clfr = RandomForestRegressor(random_state=0, bootstrap=False)
                clfr.fit(X_train, y_train)
                crsvl.append(clfr.score(X_test, y_test))
                    
        sv[' '.join(single_cat)]=np.mean(crsvl)
        sv[' '.join(single_cat)+'_std']=np.std(crsvl)
        #group and out
        #Convert dict to dataframe and choose colors
        scores = pd.DataFrame(list(sv.items()))
        scores=scores.set_index(scores[0])
        scores = scores.drop([0], 1)
        scores.columns = ['Scores']
        scores.sort_values(['Scores'], ascending = [False], inplace = True)
        scores_mean=scores.T[[x for x in scores.index if 'std' not in x]].T
        scores_mean_std=scores.T[[x for x in scores.index if 'std' in x]].T
        scores_mean_std.index=[x.replace('_std','') for x in scores_mean_std.index]
        scores=pd.concat([scores_mean,scores_mean_std],axis=1)
        scores.columns=['Mean Matrix Completion (RF)','std Matrix Completion (RF)']
        return scores,otulearn

    else:
        #run ml
        for reduced_this in list(mapdf.columns):
            
            #mapping
            tmp_otutabledf,tmpmapdfch=match(otulearn.T,mapdf[~mapdf[reduced_this].isin(filterlist)])
            tmp_otutabledf=tmp_otutabledf.T
            mapdfch=tmpmapdfch.copy()
            tmpmapdf_encode=encode_mapping(tmpmapdfch)
            
            crsvl=[] #store them
            sss = ShuffleSplit(n_splits=mean_count, test_size=test_split, random_state=0)
            sss.get_n_splits(tmp_otutabledf.as_matrix().copy().T, np.array(tmpmapdf_encode.as_matrix()))
            for train_index, test_index in sss.split(tmp_otutabledf.as_matrix().copy().T, np.array(tmpmapdf_encode.as_matrix())):
                X_train, X_test = tmp_otutabledf.as_matrix().copy().T[train_index], tmp_otutabledf.as_matrix().copy().T[test_index]
                y_train, y_test = np.array(tmpmapdf_encode.as_matrix())[train_index].T[split_tmp], np.array(tmpmapdf_encode.as_matrix())[test_index].T[split_tmp]


                if len(set(mapdf[reduced_this]))<=1: # can not learn classifiers with one label
                    break

                # learn 
                if all(isinstance(item, str) for item in mapdf[reduced_this]) or all(isinstance(item, bool) for item in mapdf[reduced_this]):
                    clfr = RandomForestClassifier(random_state=0, bootstrap=False)
                    clfr.fit(X_train, y_train)
                    crsvl.append(clfr.score(X_test, y_test))
                else: 
                    clfr = RandomForestRegressor(random_state=0, bootstrap=False)
                    clfr.fit(X_train, y_train)
                    crsvl.append(clfr.score(X_test, y_test)) 
                
            if not crsvl:
                continue
            else:
                sv[reduced_this]=np.mean(crsvl)
                sv[reduced_this+'_std']=np.std(crsvl)
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
        if 'std' not in get_s:
            new_names.append(get_s+' (n='+str(len(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))+')'+' (labels='+str(len(list(set(mapdf[~mapdf[get_s].isin(filterlist)][get_s]))))+')')
        else:
            new_names.append(get_s[:-4]+'_std'+' (n='+str(len(mapdf[~mapdf[get_s[:-4]].isin(filterlist)][get_s[:-4]]))+')'+' (labels='+str(len(list(set(mapdf[~mapdf[get_s[:-4]].isin(filterlist)][get_s[:-4]]))))+')')
    scores.index=new_names
    #group and out
    scores_mean=scores.T[[x for x in scores.index if 'std' not in x]].T
    scores_mean_std=scores.T[[x for x in scores.index if 'std' in x]].T
    scores_mean_std.index=[x.replace('_std','') for x in scores_mean_std.index]
    scores=pd.concat([scores_mean,scores_mean_std],axis=1)
    scores.columns=['Mean Matrix Completion (RF)','std Matrix Completion (RF)']

    return scores,otulearn


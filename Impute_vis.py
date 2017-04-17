import brewer2mpl
import matplotlib
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from gneiss.util import match

class PCA_niche(object):

    '''''
    This contains visuals for DEICODE
    '''''

    def niche_visual(otu,sfi,tax_index,bact_to_show,niche_plot,Pc_plot,mappingdf):
        

        #imputed
        observed_table_sfi=pd.DataFrame(sfi,index=otu.index.values,columns=otu.T.index.values) #re-order to origonal otu table
        observed_table_sfi['new_index']=tax_index # add taxa column
        otu_taxa=observed_table_sfi.set_index('new_index') # reindex
        plot_var_otu=otu_taxa.copy() # copy with taxa names
        plot_var_otu=plot_var_otu.groupby(plot_var_otu.index).sum() #group by name
        plot_var_otu=plot_var_otu.T #transpose for PCA
        
        #org
        otuch=otu.copy()
        otuch['new_index']=tax_index # add taxa column
        otuch_taxa=otuch.set_index('new_index') # reindex
        plot_var_otuch=otuch_taxa.copy() # copy with taxa names
        plot_var_otuch=plot_var_otuch.groupby(plot_var_otuch.index).sum() #group by name
        plot_var_otuch=plot_var_otuch.T #transpose for PCA
        
        pca_model=PCA(n_components=3) #WPCA
        X2=plot_var_otu.as_matrix() #pd to matrix
        X_reduced_var = pca_model.fit_transform(X2) #transform
        pccompdf = pd.DataFrame(pca_model.components_,columns=plot_var_otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
        pccompdf.sort_values([Pc_plot], ascending = [False], inplace = True) #sort weights
        highest_var_bact=list(pccompdf.index.values) # re order table by wieghts
        
        plot_var_otuch=plot_var_otuch.T
        plot_var_otuch=plot_var_otuch.reindex(highest_var_bact)
        plot_var_otuch=plot_var_otuch.T
        plot_var_otu, mapping_plot_var_otu = match(plot_var_otuch, mappingdf)
        mapping_plot_var_otu["IDs"]=list(mapping_plot_var_otu.index.values)
        
        # add sorted ids into dict by desired metadata factor
        nichesetdict=mapping_plot_var_otu.set_index('IDs')[niche_plot].to_dict()
        nichesetdict=dict([a, x] for a, x in nichesetdict.items())
        # sort ID's into niche
        grouped_nichesetdict = {}
        for sample_ids, group_value in sorted(nichesetdict.items()):
            grouped_nichesetdict.setdefault(group_value, []).append(sample_ids)
        
        nich_col=[]
        for check_samp_niche in list(plot_var_otu.index.values):
            for niche_group, all_sample_id_check in grouped_nichesetdict.items():
                for sample_id_check in all_sample_id_check:
                    if check_samp_niche==sample_id_check:
                        nich_col.append(niche_group)
        
        plot_var_otu["niche"]=nich_col
        plot_var_otu_grouped = plot_var_otu.groupby('niche').describe()
        out_niche_linked = {}
        index_rename=[]
        index_mean=[]
        index_std=[]
        for niche in list(grouped_nichesetdict.keys()):
            for bact_var in highest_var_bact[:bact_to_show]:
                out_niche_linked.setdefault(niche, []).append(plot_var_otu_grouped[bact_var][niche]['max'])
                out_niche_linked.setdefault(niche, []).append(plot_var_otu_grouped[bact_var][niche]['std'])
        for bact_var in highest_var_bact[:bact_to_show]:
            index_rename.append('%s'%bact_var)
            index_rename.append('%s_std'%bact_var)
            index_mean.append('%s'%bact_var)
            index_std.append('%s_std'%bact_var)
    
        out_niche_linkeddf=pd.DataFrame(out_niche_linked, index=index_rename)

        return out_niche_linkeddf,observed_table_sfi,index_mean,index_std,highest_var_bact,pccompdf


    def plot_niche(out_niche_linkeddf,observed_table_sfi,mappingdf,encoded_mapping,niche_plot,Pc_plot,index_mean,index_std,le,cont):

        #PCA
        observed_table_sfi, mapping = match(observed_table_sfi.T, mappingdf)
        pca_model=PCA(n_components=2)
        X2=observed_table_sfi.as_matrix()
        X_reduced2 = pca_model.fit_transform(X2)
        
    
        if Pc_plot=="PC-1":
        
            if cont==True:
                
                index_mean_new=[]
                index_std_new=[]
                new_taxa=[]
                #lowest taxaclassification + (phylum)
                for nm in out_niche_linkeddf.index.values.tolist():
                    #find lowest taxa level
                    q=0
                    while q<(len(nm.split(';'))-1):
                        lw=nm.split(';')[q]
                        q+=1
                        if len(nm.split(';')[q])==3:
                            break
                    if '[' in lw.split('__')[1]:
                        lw=lw.split('__')[1]
                        lw=lw[1:-1]
                    else:
                        lw=lw.split('__')[1]

                    if 'std' in nm:
                        index_std_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                    else:
                        index_mean_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                out_niche_linkeddf.index=new_taxa
                #plot
                fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(20, 18),sharey=False)
                out_niche_linkeddf.T[index_mean_new[:-5]].plot(kind='area',lw=6, alpha = 0.92,rot=0,colormap="Set1",sharey=True,fontsize=15,ax=ax1)
                ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
                Y=encoded_mapping[niche_plot][0].tolist()
                ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(list(mappingdf[niche_plot])),cmap=plt.cm.cool,s=200)
                p=ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(list(mappingdf[niche_plot])),cmap=plt.cm.cool,s=200)
                fig.colorbar(p,orientation='horizontal')
                # Set common labels
                ax1.set_xlabel(('Mean Frequency in %s Niche'%niche_plot))
                ax1.set_ylabel('Niche')
                ax2.set_xlabel('PC-1')
                ax2.set_ylabel('PC-2')
                ax1.set_title(('OTUs With Highest Weighted Varaince in %s'%niche_plot))
                ax2.set_title(('PCA on Fully Dense Imputed OTU Table (colored by %s)'%niche_plot))
                fig.set_tight_layout(True)
                return plt

            if cont==False:
                
                ids=['PC1','PC2']
                spPCAplot=pd.DataFrame(X_reduced2,observed_table_sfi.index.values,ids)
                spPCAplot['labels']=mappingdf[niche_plot]
                #clean
                index_mean_new=[]
                index_std_new=[]
                new_taxa=[]
                #lowest taxaclassification + (phylum)
                for nm in out_niche_linkeddf.index.values.tolist():
                    #find lowest taxa level
                    q=0
                    while q<(len(nm.split(';'))-1):
                        lw=nm.split(';')[q]
                        q+=1
                        if len(nm.split(';')[q])==3:
                            break
                    if '[' in lw.split('__')[1]:
                        lw=lw.split('__')[1]
                        lw=lw[1:-1]
                    else:
                        lw=lw.split('__')[1]
                    
                    if 'std' in nm:
                        index_std_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                    else:
                        index_mean_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                out_niche_linkeddf.index=new_taxa
                #plot
                fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(20, 18),sharey=False)
                #scatter
                sns.swarmplot(x="PC1", y="PC2", data=spPCAplot, hue="labels", size=10,ax=ax1)
                ax1.set_xlabel('$PC-1$')
                ax1.set_ylabel('$PC-2$')
                ax1.legend(loc=2,prop={'size':22},bbox_to_anchor=(.59, 1.03))
                #barh
                out_niche_linkeddf.T[index_mean_new].plot(kind='barh',width=.9, xerr=out_niche_linkeddf.T[index_std_new].values.T, alpha = 0.92,error_kw=dict(elinewidth=1,capsize=2,barsabove=True,ecolor='k',ms=1, mew=1),rot=0,colormap="Set1",sharey=True,fontsize=15,ax=ax2)
                ax2.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
                fig.set_tight_layout(True)
                return plt

        if Pc_plot=="PC-2":


            if cont==True:
    
                index_mean_new=[]
                index_std_new=[]
                new_taxa=[]
                #lowest taxaclassification + (phylum)
                for nm in out_niche_linkeddf.index.values.tolist():
                    #find lowest taxa level
                    q=0
                    while q<(len(nm.split(';'))-1):
                        lw=nm.split(';')[q]
                        q+=1
                        if len(nm.split(';')[q])==3:
                            break
                    if '[' in lw.split('__')[1]:
                        lw=lw.split('__')[1]
                        lw=lw[1:-1]
                    else:
                        lw=lw.split('__')[1]
                    
                    if 'std' in nm:
                        index_std_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                    else:
                        index_mean_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                out_niche_linkeddf.index=new_taxa
                #plot
                fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(20, 18),sharey=False)
                out_niche_linkeddf.T[index_mean_new[:-5]].plot(kind='area',lw=6, alpha = 0.92,rot=0,colormap="Set1",sharey=True,fontsize=15,ax=ax1)
                ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
                Y=encoded_mapping[niche_plot][0].tolist()
                ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(list(mappingdf[niche_plot])),cmap=plt.cm.cool,s=200)
                p=ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(list(mappingdf[niche_plot])),cmap=plt.cm.cool,s=200)
                fig.colorbar(p,orientation='horizontal')
                # Set common labels
                ax1.set_xlabel(('Mean Frequency in %s Niche'%niche_plot))
                ax1.set_ylabel('Niche')
                ax2.set_xlabel('PC-1')
                ax2.set_ylabel('PC-2')
                ax1.set_title(('OTUs With Highest Weighted Varaince in %s'%niche_plot))
                ax2.set_title(('PCA on Fully Dense Imputed OTU Table (colored by %s)'%niche_plot))
                fig.set_tight_layout(True)
                return plt

            if cont==False:
                
                ids=['PC1','PC2']
                spPCAplot=pd.DataFrame(X_reduced2,observed_table_sfi.index.values,ids)
                spPCAplot['labels']=mappingdf[niche_plot]
                #clean
                index_mean_new=[]
                index_std_new=[]
                new_taxa=[]
                #lowest taxaclassification + (phylum)
                for nm in out_niche_linkeddf.index.values.tolist():
                    #find lowest taxa level
                    q=0
                    while q<(len(nm.split(';'))-1):
                        lw=nm.split(';')[q]
                        q+=1
                        if len(nm.split(';')[q])==3:
                            break
                    if '[' in lw.split('__')[1]:
                        lw=lw.split('__')[1]
                        lw=lw[1:-1]
                    else:
                        lw=lw.split('__')[1]

                    if 'std' in nm:
                        index_std_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')'+'std')
                    else:
                        index_mean_new.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                        new_taxa.append(lw+' ('+nm.split(';')[0].split('__')[1]+')')
                out_niche_linkeddf.index=new_taxa
                #plot
                fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(13, 8),sharey=False)
                #scatter
                sns.swarmplot(x="PC1", y="PC2", data=spPCAplot, hue="labels", size=10,ax=ax1)
                ax1.set_xlabel('$PC-1$')
                ax1.set_ylabel('$PC-2$')
                ax1.legend(loc=2,prop={'size':22},bbox_to_anchor=(.59, 1.03))
                #barh
                out_niche_linkeddf.T[index_mean_new].plot(kind='barh',width=.9, xerr=out_niche_linkeddf.T[index_std_new].values.T, alpha = 0.92,error_kw=dict(elinewidth=1,capsize=2,barsabove=True,ecolor='k',ms=1, mew=1),rot=0,colormap="Set1",sharey=True,fontsize=15,ax=ax2)
                ax2.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
                fig.set_tight_layout(True)
                return plt

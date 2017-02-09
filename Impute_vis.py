import brewer2mpl
import matplotlib
import pandas as pd
from sklearn.decomposition import PCA
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
        plot_var_otu=plot_var_otu.groupby(plot_var_otu.index).mean() #group by name
        plot_var_otu=plot_var_otu.T #transpose for PCA
        
        #org
        otuch=otu.copy()
        otuch['new_index']=tax_index # add taxa column
        otuch_taxa=otuch.set_index('new_index') # reindex
        plot_var_otuch=otuch_taxa.copy() # copy with taxa names
        plot_var_otuch=plot_var_otuch.groupby(plot_var_otuch.index).mean() #group by name
        plot_var_otuch=plot_var_otuch.T #transpose for PCA
        
        pca_model=PCA(n_components=3) #PCA
        X2=plot_var_otu.as_matrix() #pd to matrix
        X_reduced_var = pca_model.fit_transform(X2) #transform
        pccompdf = pd.DataFrame(pca_model.components_,columns=plot_var_otu.columns,index = ['PC-1','PC-2','PC-3']).T #get wieghts
        pccompdf.sort_values([Pc_plot], ascending = [False], inplace = True) #sort weights
        highest_var_bact=list(pccompdf.index.values) # re order table by wieghts
        
        plot_var_otuch=plot_var_otuch.T
        plot_var_otuch['new_index']=highest_var_bact
        plot_var_otuch = plot_var_otuch.set_index('new_index')
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
                out_niche_linked.setdefault(niche, []).append(plot_var_otu_grouped[bact_var][niche]['mean'])
                out_niche_linked.setdefault(niche, []).append(plot_var_otu_grouped[bact_var][niche]['std'])
        for bact_var in highest_var_bact[:bact_to_show]:
            index_rename.append('%s_mean'%bact_var)
            index_rename.append('%s_std'%bact_var)
            index_mean.append('%s_mean'%bact_var)
            index_std.append('%s_std'%bact_var)
    
        out_niche_linkeddf=pd.DataFrame(out_niche_linked, index=index_rename)

        return out_niche_linkeddf,observed_table_sfi,index_mean,index_std,highest_var_bact,pccompdf


    def plot_niche(out_niche_linkeddf,observed_table_sfi,mappingdf,encoded_mapping,niche_plot,Pc_plot,index_mean,index_std,le,cont):

        #PCA
        observed_table_sfi, mapping = match(observed_table_sfi.T, mappingdf)
        pca_model=PCA(n_components=3)
        X2=observed_table_sfi.as_matrix()
        X_reduced2 = pca_model.fit_transform(X2)
        
        if cont==True:
            
            fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(25, 10),sharey=False)
            out_niche_linkeddf.T[index_mean].plot(kind='bar',width=.78, yerr=out_niche_linkeddf.T[index_std].values.T, alpha = 0.92,error_kw=dict(elinewidth=1,capsize=2,barsabove=True,ecolor='k',ms=1, mew=1),rot=0,colormap="Set3",sharey=True,fontsize=15,ax=ax1)
            ax1.legend(loc=2,prop={'size':16},bbox_to_anchor=(1.0, 1.0))
            Y=encoded_mapping[niche_plot][0].tolist()
            ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
            p=ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=list(encoded_mapping[bestclassifier][0]),cmap=plt.cm.cool,s=200)
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
            
            fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(25, 25),sharey=False)
            
            #bar
            out_niche_linkeddf.T[index_mean].plot(kind='bar',width=.78, yerr=out_niche_linkeddf.T[index_std].values.T, alpha = 0.92,error_kw=dict(elinewidth=1,capsize=2,barsabove=True,ecolor='k',ms=1, mew=1),rot=0,colormap="Set3",sharey=True,fontsize=30,ax=ax1)
            ax1.legend(loc=2,prop={'size':25},bbox_to_anchor=(1.0, 1.0))
            Y=mapping[niche_plot].tolist()
            
            #scatter (with legend for non coninous data)
            Y=list(map(int, encoded_mapping[niche_plot][0]))
            Y2=list(encoded_mapping[niche_plot][1])
            bmap7 = brewer2mpl.get_map('Set1','qualitative',9,reverse=True)
            colors = bmap7.mpl_colors
            handles, labels = ax2.get_legend_handles_labels()
            ax2.legend(handles, labels)
            scatter_proxy=[]
            k=set(Y2)
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
                if key==int((len(set(Y2)))):
                    for q in value:
                        scatter_proxy.append(matplotlib.lines.Line2D([0],[0], linestyle="none", c=colors[q],markersize=30, marker = 'o'))
            
            if int(len(k))>=3:
                numrows=4
            else:
                numrows=int(len(k))
            ax2.legend(scatter_proxy, k, numpoints = 1,loc=2,prop={'size':40},bbox_to_anchor=(1.0, 1.0),fontsize = 'large', labelspacing=.1, borderaxespad=0.)
            ax2.scatter(X_reduced2[:, 0], X_reduced2[:, 1], c=Y,cmap=plt.cm.Set1_r,s=350)
            # Set common labels
            ax1.set_xlabel(('Mean Frequency in %s Niche'%niche_plot), fontsize=20)
            ax1.set_ylabel('Niche (mean frequency)', fontsize=20)
            ax2.set_xlabel('PC-1', fontsize=20)
            ax2.set_ylabel('PC-2', fontsize=20)
            ax1.set_title(('OTUs With Highest Weighted Varaince in %s'%niche_plot), fontsize=35)
            ax2.set_title(('PCA on Fully Dense Imputed OTU Table (colored by %s)'%niche_plot), fontsize=35)
            fig.set_tight_layout(True)
            return plt

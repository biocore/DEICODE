#utils
from __future__ import division
import pandas as pd
import numpy as np
from biom import load_table
import sys


def matchtable(df1,df2):

    '''
    match dataframes 
    '''
    #check type
    df1.columns=map(str, list(df1.columns))
    df2.columns=map(str, list(df2.columns))
    #match
    df1=df1[list(set(list(df1.columns)) & set(list(df2.columns)))]
    df2=df2[list(set(list(df1.columns)) & set(list(df2.columns)))]
    return df1,df2

def clean_names_match_map(df,dfmap,num=13):
    
    '''
    input: cleans names for otu table and matches it to mapping
    
    output: mtached mapping and otu tables
    
    '''
    
    df.columns=[x[num:].split('DNA')[0] for x in list(df.columns)]
    return matchtable(df,dfmap)
        
def convert_biom_to_pandas(table): # covert biom
    
    feature_table = pd.DataFrame(np.array(table.matrix_data.todense()).T,index=table.ids(axis='sample'),columns=table.ids(axis='observation'))
    feature_ids = table.ids(axis='observation')
    mapping = {i: table.metadata(id=i, axis='observation')['taxonomy'] for i in feature_ids}
    for key, value in mapping.items():
        nvalue = ';'.join(value[1:])
        mapping.update({key:nvalue})
    taxonomy = pd.DataFrame(mapping, index=['taxonomy']).T
    
    return feature_table, taxonomy

def taxa_lvl(taxa_names,txlvl):
    
    return [";".join(t.split(";")[:txlvl]) for t in taxa_names]

def dfimport(in_biom,map_file,filter_count=0,txlvl=6):
    
    
    '''
    input: path to otu table (can be .biom or tabbed) and mapping 
    
    output: matched dataframes for mapping and otu table and taxonomy names for level specified 
    
    '''
    #mapping file tab delim
    mappingdf= pd.read_table('%s'%map_file, index_col=0,low_memory=False)
    mappingdf=mappingdf.replace(np.nan,'Unknown', regex=True)

    #otu table (biom or csv)
    try:
        filename=in_biom.split('/')[-1]
    except:
        filename=in_biom
    if filename.split('.')[-1]=="biom":
        #BIOM
        #load table
        table = load_table('%s'%in_biom)
        read_filter = lambda val, id_, md: sum(val) > filter_count
        table.filter(read_filter, axis='sample')
        table.filter(read_filter, axis='observation')
        otu, taxanames = convert_biom_to_pandas(table)
        otu=otu.T
        otu.fillna(0)
        taxanames.index=[str("OTU_%s"%str(q)) for q in range(0,len(otu.index.values))]
        if min(otu.shape)<=1:
            sys.exit('Import error less than two samples or features: please check that your data is tab delimited or in biom file format')
        #tax to level (i.e. 6 is usually species)
        otu.index=[str("OTU_%s"%str(q)) for q in range(0,len(otu.index.values))]
        otu,mappingdf=matchtable(otu,mappingdf.T)
        return otu,mappingdf,taxanames
    elif filename.split('.')[-1]=="csv" or filename.split('.')[-1]=="tsv" or filename.split('.')[-1]=="txt":
        #csv
        otu=pd.read_table('%s'%in_biom, index_col=0)
        otu.fillna(0)
        taxanames = pd.DataFrame(taxa_lvl(list(otu.index),txlvl), index=[str("OTU_%s"%str(q)) for q in range(0,len(otu.index.values))]).T
        if min(otu.shape)<=1:
            sys.exit('Import error less than two samples or features: please check that your data is tab delimited or in biom file format')
        #tax to level (i.e. 6 is usually species)
        otu.index=[str("OTU_%s"%str(q)) for q in range(0,len(otu.index.values))]
        otu,mappingdf=matchtable(otu,mappingdf.T)
        return otu,mappingdf,taxanames
    else:
        sys.exit('Import error: please check that your data is one of the following file formats (.csv,.biom,.txt,.tsv)')



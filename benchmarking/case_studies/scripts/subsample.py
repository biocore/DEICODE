import pandas as pd
import numpy as np
import os
from gneiss.util import match
from biom.util import biom_open
from biom import load_table
from collections import Counter
import numpy as np
import shutil
import os
np.random.seed(42)

def make_safe_dir(tmp_w):
    if not os.path.exists(tmp_w):
        os.makedirs(tmp_w)
    else:
        shutil.rmtree(tmp_w)           
        os.makedirs(tmp_w)

# empty dict to save each case study in 
case_study={}
case_study['Sponges']={}
case_study['Sleep_Apnea']={}
case_study['Sponges']['factor']='health_status'
case_study['Sleep_Apnea']['factor']='exposure_type'

for dataset_ in [dir_ for dir_ in os.listdir('data/') if dir_!='.DS_Store']:
    #import and filter samples for >1000 reads/sample
    in_biom = os.path.join('data',dataset_,'table.biom')
    table = load_table(in_biom)
    sample_filter = lambda val, id_, md: sum(val) > 1000
    table = table.filter(sample_filter, axis='sample')
    # Get metadata for table matched
    mapping_df = pd.read_table(os.path.join('data',dataset_,'metadata.txt')
                            ,index_col=0, low_memory=False)
    keep_ = list(set(table.ids())&set(mapping_df.index))
    # save tables
    case_study[dataset_]['biom_table'] = table.copy().filter(keep_,axis='sample')
    case_study[dataset_]['metadata'] = mapping_df.reindex(index=keep_)
# sub-sample each dataset
for dataset_ in case_study.keys():
    sub_sample={}
    sub_sample_table={} 
    # repeat randomized sub-sample 10 fold times
    for fold_ in range(1,10):
        # Randomly Sub-sample to balance groups
        balance_tmp=case_study[dataset_]['metadata'].copy()
        balance_tmp_table=case_study[dataset_]['biom_table'].copy()
        # no-subsample find count and save
        count=max(list(Counter(balance_tmp[case_study[dataset_]['factor']]).values()))
        sub_sample[(fold_,int(count*2))]=balance_tmp
        #save table 
        sub_sample_table[(fold_,int(count*2))] = balance_tmp_table.copy().filter(balance_tmp.index, 
                                                                          axis='sample')

        balance_tmp['org_index']=balance_tmp.index
        balance_tmp=balance_tmp.groupby(case_study[dataset_]['factor'])
        balance_tmp=balance_tmp.apply(lambda x: x.sample(balance_tmp.size().min()).reset_index(drop=True))
        balance_tmp.index=balance_tmp.index.droplevel()    
        start=list(Counter(balance_tmp[case_study[dataset_]['factor']]).values())[0]
        # save balanced groups
        balance_tmp.index=balance_tmp['org_index']
        sub_sample[(fold_,int(start*2))]=balance_tmp
        #save table 
        sub_sample_table[(fold_,int(start*2))] = balance_tmp_table.copy().filter(balance_tmp.index, 
                                                                          axis='sample')

        sub_group=list(set(balance_tmp[case_study[dataset_]['factor']]))[0]
        per_group_sub=np.arange(15,start,20)
        idx = np.round(np.linspace(0, len(per_group_sub) - 1, 5)).astype(int)
        per_group_sub=per_group_sub[idx]
        balance_tmp_re=balance_tmp.copy()
        # sub-sample reads randomly for on group then re-balance
        for sub_count in per_group_sub[::-1]:
            # get sub-index for sub-group
            sub_index = balance_tmp[balance_tmp[case_study[dataset_]['factor']].isin([sub_group])]
            drop_indices = np.random.choice(sub_index.index, 
                                            start-sub_count, 
                                            replace=False)
            balance_tmp = balance_tmp.drop(drop_indices)
            # re-balance after removing (start-sub_count) indecies randomly
            balance_tmp = balance_tmp.groupby(case_study[dataset_]['factor'])
            balance_tmp = balance_tmp.apply(lambda x: x.sample(balance_tmp.size().min()).reset_index(drop=True))
            balance_tmp.index = balance_tmp.index.droplevel()    
            start = list(Counter(balance_tmp[case_study[dataset_]['factor']]).values())[0]
            # save balanced groups
            balance_tmp.index = balance_tmp['org_index']
            sub_sample[(fold_,int(sub_count*2))] = balance_tmp
            #save table 
            sub_sample_table[(fold_,int(sub_count*2))] = balance_tmp_table.copy().filter(balance_tmp.index, 
                                                                              axis='sample')

    for (fold_out,sample_N_),metadata_tmp in sub_sample.items():
        dir_tmp_='sub_sample/biom_tables_'+str(dataset_)+'/'+str(fold_out)+'_'+str(sample_N_)
        make_safe_dir(dir_tmp_)
        metadata_tmp.to_csv(dir_tmp_+'/metadata.tsv',sep='\t')

    for (fold_out,sample_N_),out_table_tmp in sub_sample_table.items():
        dir_tmp_='sub_sample/biom_tables_'+str(dataset_)+'/'+str(fold_out)+'_'+str(sample_N_)
        name_=dir_tmp_+'/table.biom'
        with biom_open(name_, 'w') as write_table:
            out_table_tmp.remove_empty().to_hdf5(write_table,'subsample_'+str(dataset_))

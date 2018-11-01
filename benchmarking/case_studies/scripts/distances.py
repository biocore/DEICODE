import pandas as pd
import numpy as np
import os
from gneiss.util import match
from biom.util import biom_open
from biom import load_table
from collections import Counter
import numpy as np
import subprocess
import shutil 

def make_safe_dir(tmp_w):
    if not os.path.exists(tmp_w):
        os.makedirs(tmp_w)
    else:
        shutil.rmtree(tmp_w)           
        os.makedirs(tmp_w)

# empty dict to save each case study in 
dir_path = os.getcwd()
#ranks
ranks={}
ranks['Sleep_Apnea'] = 3 # three clusters here
ranks['Sponges'] = 2 # two clusters here 
for dataset_ in [dir_ for dir_ in os.listdir('data/') if dir_!='.DS_Store']:
    tree_path = os.path.join(dir_path,'data',dataset_,'tree_relabelled.tre')
    tree_path_qiime = os.path.join(dir_path,'data',dataset_,'tree_relabelled.qza')
    # make qiime2 tree
    subprocess.call("qiime tools import --input-path %s --type 'Phylogeny[Rooted]' --output-path %s"%(tree_path,tree_path_qiime), shell=True)
    subpath_=os.path.join('sub_sample','biom_tables_'+dataset_)
    for sub_set in [dir_ for dir_ in os.listdir(subpath_) if dir_!='.DS_Store']:
        sub_path_sub = os.path.join(dir_path,subpath_,sub_set+'/')
        table_path = os.path.join(dir_path,subpath_,sub_set,'table.biom')
        metadata_path_qiime = os.path.join(dir_path,subpath_,sub_set,'metadata.tsv')  
        # make qiime2 table form biom 
        table_path_qiime=os.path.join(dir_path,subpath_,sub_set,'table.qza')
        subprocess.call("qiime tools import --input-path %s --type 'FeatureTable[Frequency]' --source-format BIOMV210Format --output-path %s"%(table_path,table_path_qiime), shell=True)
        # sequencing depth cleaning (1000 read/sample)
        subprocess.call("qiime feature-table filter-samples --i-table %s --p-min-frequency 1000 --o-filtered-table %s"%(table_path_qiime,table_path_qiime), shell=True)
        # trim biom of any tree issues (should not be any but why not)
        subprocess.call("qiime phylogeny filter-table --i-table %s --i-tree %s --o-filtered-table %s"%(table_path_qiime,tree_path_qiime,table_path_qiime), shell=True)
        # rarefy (1,000)
        subprocess.call("qiime feature-table rarefy --i-table %s --p-sampling-depth 1000 --o-rarefied-table %s"%(table_path_qiime,table_path_qiime), shell=True)
        # bray-curtis
        bray_dist = os.path.join(dir_path,subpath_,sub_set,'Bray_Distance.qza')
        subprocess.call("qiime diversity beta --i-table %s --p-metric 'braycurtis' --o-distance-matrix %s"%(table_path_qiime,bray_dist), shell=True)
        # generalized weighted alpha=1.0
        wone_dist = os.path.join(dir_path,subpath_,sub_set,'GUniFrac_alpha_one_Distance.qza')
        subprocess.call("qiime diversity beta-phylogenetic-alt --i-table %s --i-phylogeny %s --p-metric generalized_unifrac --p-alpha 1.0 --o-distance-matrix %s"%(table_path_qiime,tree_path_qiime,wone_dist), shell=True)           
        # extract txt files from each distance
        for convert_ in [bray_dist,wone_dist]:
            subprocess.call("qiime tools export %s --output-dir %s"%(convert_,convert_.replace('.qza','')),shell=True)      
        for convert_ in [bray_dist,wone_dist]:
            start_=os.path.join(convert_.replace('.qza',''),'distance-matrix.tsv')
            end_=convert_.replace('.qza','.tsv')
            shutil.move(start_,end_)
            os.rmdir(convert_.replace('.qza',''))
        # remove temp qza files to save space
        rm_tmp_files=[os.path.join(sub_path_sub,file_) for file_ in os.listdir(sub_path_sub) if '.qza' in file_]
        _=[os.remove(file_) for file_ in rm_tmp_files] 
        # run RPCA
        subprocess.call("deicode_rpca --in_biom %s --output_dir %s --rank %s --min_sample_depth 1000"%(table_path,sub_path_sub,str(ranks[dataset_])), shell=True)                

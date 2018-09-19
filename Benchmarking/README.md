# Simulations

The results and visualizations can be found in the jupyter notebook. 

### Building simulations on real data 
 
* Clusters - The Keyboard Dataset (Qiita ID 232) - cluster_build.py


We will use a study of keyboard skin communities done by Noah Fierer (University of Colorado) as a model dataset; they discover that

    "diversity of skin-associated bacterial communities is far higher than previously recognized, with a high degree of inter-individual variability in the composition of bacterial communities." 

The study goes on describing the amazing forensic implications of matching a person's skin microbial community to a recently interacted surface. However, in the process a prime example of a low-rank microbial dataset was created and has recently become a favorite benchmarking dataset for new tools. Here we will fit a model on this data and generate data with different numbers of cluster (ranks) and separation (by overlap between clusters).


* Gradients - The 88 Soils Dataset (Qiita ID 103) - gradient_build.py   

We will use the 88 soils dataset as a model dataset for gradients. 88 soils dataset is a study done by Noah Fierer (University of Colorado) they discover that

    "We found that overall bacterial community composition, as measured by pairwise UniFrac distances, was significantly correlated with differences in soil pH (r = 0.79), largely driven by changes in the relative abundances of Acidobacteria, Actinobacteria, and Bacteroidetes across the range of soil pHs. In addition, soil pH explains a significant portion of the variability associated with observed changes in the phylogenetic structure within each dominant lineage. The overall phylogenetic diversity of the bacterial communities was also correlated with soil pH (R2 = 0.50), with peak diversity in soils with near-neutral pHs. Together, these results suggest that the structure of soil bacterial communities is predictable, to some degree, across larger spatial scales, and the effect of soil pH on bacterial community composition is evident at even relatively coarse levels of taxonomic resolution."

Here we will explain and investigate the latent gradient structure of the microbial communities between the different samples. Then we will fit a model on this data and generate data with different band widths (sigma). The output of gradient_build.py is a csv file with simulation with added noise and sparsity along with a fully dense base truth dataset. 


### Compare methods with KL-Divergence on clr-transformed data

From the modeled output from cluster_build.py and gradient_build.py we will fit multiple completion methods. These include:

1. Pseudo Counts (one count applied to all data)
2. Soft Impute 
3. Iterative SVD
4. KNN
5. OptSpace

The output is then compared to the base truth data by the KL-Divergence of the clr transform. This will be done through cluster_results.py & gradient_results.py the output will be a result.csv file in the respective directories.  


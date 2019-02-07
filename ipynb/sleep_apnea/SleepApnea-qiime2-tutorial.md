# Tutorial
 
**Note**: This guide assumes you have installed QIIME 2 using one of the procedures in the [install documents](https://docs.qiime2.org/2019.1/install/) and have installed [q2-deicode](https://library.qiime2.org/plugins/q2-deicode).

## Introduction 

In this tutorial you will learn how to interpret and perform Robust Aitchison PCA through QIIME 2. The focus of this tutorial is compositional beta diversity. There are many beta diversity metrics that have been proposed, all with varying benefits on varying data structures. However, presence/absence metric often prove to give better results than those that rely on abundances (i.e. unweighted vs. weighted UniFrac). One component of this phenomenon is that the interpretation of relative abundances can provide spurious results (see [the differential abundance analysis introduction](https://docs.qiime2.org/2019.1/tutorials/gneiss/)). One solution to this problem is to use a compositional distance metric such as [Aitchison distance](https://en.wikipedia.org/wiki/Aitchison_geometry). 


As a toy example let’s build three taxa. These three taxa represent common distributions we see in microbiome datasets. Where the first taxon is increasing exponentially across samples, this is a trend that we would be interested in. However, taxon 2 and 3 have much higher counts and taxon 3 is randomly fluctuating across samples.  

![](https://cdn1.imggmi.com/uploads/2019/2/5/ed98350fbb7df22df074c3751268ad09-full.png)

In our distances below we have Euclidean, Bray-Curtis, Jaccard, and Aitchison distances (from left to right). We can see that the abundance based metrics Euclidean and Bray-Curtis are heavily influenced by the abundance of taxon 3 and seem to randomly fluctuate. In the presence/absence metric, Jaccard, we see that the distance saturates to one very quickly. However, in the Aitchison distance we see a linear curve representing taxon 1. The reason the distance is linear is because Aitchison distance relies on log transforms (the log of the exponential trend of taxon 1 is linear). 


![](https://cdn1.imggmi.com/uploads/2019/2/5/ccf5feb1e1cfeda1329689abe949a3c7-full.png)

From this toy example, it is clear that Aitchison distance better accounts for the proportions. However, we made the unrealistic assumption in our toy example that there were no zero counts. In real microbiome datasets there are a large number of zeros (i.e. sparsity). Sparsity complicates log ratio transformations because the log-ratio of zero is undefined. To solve this, pseudo counts are often used but that can often skew results (see [Naught all zeros in sequence count data are the same](https://www.biorxiv.org/content/10.1101/477794v1)). 

Robust Aitchison PCA solves this problem in two steps:

**1.** Compostional preprocessing using the centered log ratio transform on only the non-zero values of the data (no pseudo count)

![](https://latex.codecogs.com/gif.latex?rclr%28x%29%20%3D%20%5B%5Clog%5Cfrac%7Bx_%7B1%7D%7D%7Bg_%7Br%7D%28x%29%7D%20%2C%20...%20%2C%20%5Clog%5Cfrac%7Bx_%7BD%7D%7D%7Bg_%7Br%7D%28x%29%7D%5D)

![](https://latex.codecogs.com/gif.latex?g_%7Br%7D%28x%29%20%3D%20%28%5Cprod_%7Bi%20%5Cin%20%5COmega%20_%7Bx%7D%7D%5E%7B%20%7D%20x_%7Bi%7D%29%5E%7B1/%5Cleft%20%7C%20%5COmega%20_%7Bx%7D%20%5Cright%20%7C%7D)

**2.** Dimensionality reduction through Robust PCA on only the non-zero values of the data ( [matrix completion]( https://arxiv.org/pdf/0906.2027.pdf)). 

![](https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cboldsymbol%7BU%7D%2C%20%5Cboldsymbol%7BV%7D%7D%20%5C%3B%20%5Cbigg%5Cvert%20%5CLambda%20%5Cleft%28%20%5Cboldsymbol%7BY%7D%20-%20%5Cboldsymbol%7BU%7D%20%5Cboldsymbol%7BS%7D%20%5Cboldsymbol%7BV%7D%5E%7BT%7D%20%5Cright%29%20%5Cbigg%5Cvert%20_%7B2%7D%5E%7B2%7D)

To demonstrate this in action we will run an example dataset below, where the output can be viewed as a compositional biplot through emperor. 

## Example 

In this example we will Robust Aitchison PCA via deicode on a study of sleep apnea published in [Tripathi et al. (2018)](https://msystems.asm.org/content/3/3/e00020-18). The study consists of mouse fecal samples and focuses on comparing the gut microbiome of animals exposed to intermittent hypoxia and hypercapnia (IHH; as a model of obstructive sleep apnea) to controls exposed to room air (air). 

Before beginning the tutorial let’s make a directory to store the data

```shell
mkdir qiime2-sleep-apnea-tutorial
cd qiime2-sleep-apnea-tutorial
```

#### Table
[**Download URL**](https://github.com/biocore/DEICODE/blob/master/ipynb/sleep_apnea/qiime2-sleep-apnea-tutorial/qiita_10422_table.biom.qza)
**save as:** qiita_10422_table.biom.qza 

#### Sample Metadata
[**Download URL**]( https://github.com/biocore/DEICODE/blob/master/ipynb/sleep_apnea/qiime2-sleep-apnea-tutorial/qiita_10422_metadata.tsv)
**save as:** qiita_10422_metadata.tsv

#### Feature Metadata
[**Download URL**]( https://github.com/biocore/DEICODE/blob/master/ipynb/sleep_apnea/qiime2-sleep-apnea-tutorial/taxonomy.qza)
**save as:** taxonomy.qza

Using qiita_10422_table.biom.qza, of the type raw count table (FeatureTable[Frequency]), we will generate our beta diversity ordination file. There are a few parameters to deicode that we may want to consider. The first is filtering cutoffs, these are p-min-feature-count and p-min-sample-count. Both of these parameters accept integer values and remove feature or samples, respectively, with sums below this cutoff. The feature cut-off is useful in the case that features with very low total counts among all samples represent contamination or chimeric sequences. The sample cut off is useful for the case that some sample received very few reads relative to other samples.

**Note:** it is _not_ recommended to bin your features by taxonomic assignment (i.e. by genus level). 
**Note:** it is _not_ recommended to rarefy your data before using deicode. 

The other two parameters are --p-rank and --p-iterations. These parameters should rarely have to change from the default. However, the minimum value of --p-rank can be 1 and the maximum recommended value is 10. Similarly, the minimum value of --p-iterations is 1 and is recommended to be below 500.  

Now that we understand the acceptable parameters, we are ready to run deicode.  

```shell
 qiime dev refresh-cache
```
```shell
 qiime deicode rpca \
    --i-table qiita_10422_table.biom.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot ordination.qza \
    --o-distance-matrix distance.qza
```
**Output:** PCoAResults % Properties(['biplor]) to: ordination.qza, DistanceMatrix to: distance.qza

Now that we have our ordination file, with type (PCoAResults % Properties(['biplot'])), we are ready to visualize the results. This can be done using the [emperor](https://docs.qiime2.org/2019.1/plugins/available/emperor/) biplot functionality. In this case we will include metadata for our features (optional) and our samples (required). 

```shell
qiime emperor biplot \
    --i-biplot ordination.qza \
    --m-sample-metadata-file qiita_10422_metadata.tsv \
    --m-feature-metadata-file taxonomy.qza \
    --o-visualization biplot.qzv \
    --p-number-of-features 8
```
**Output:** biplot.qzv

The interpretation of the compositional biplot may differ from classical biplot interpretation (we can view the qzv file at [view.qiime2](https://view.qiime2.org). The important features with regard to sample clusters are not a single arrow but by the log ratio between features represented by arrows pointing in different directions. A visualization tool for these log ratios is coming soon to QIIME 2. 

![](https://cdn1.imggmi.com/uploads/2019/2/6/eaf9c58ee3b00949fcd4947333376a03-full.png)

From this visualization we noticed that exposure_type seems explain the clusters well. We can run [PERMANOVA](https://docs.qiime2.org/2019.1/plugins/available/diversity/beta-group-significance/) on the distances to get a statistical significance for this. 

```shell
 qiime diversity beta-group-significance \
    --i-distance-matrix distance.qza \
    --m-metadata-file qiita_10422_metadata.tsv \
    --m-metadata-column exposure_type \
    --p-method permanova \
    --o-visualization exposure_group_significance.qzv
```

Indeed we can now see that the clusters we saw in the biplot were significant by viewing the exposure_group_significance.qzv at [view.qiime2](https://view.qiime2.org).

![](https://cdn1.imggmi.com/uploads/2019/2/6/f38d41bce26d8d7930db270680921130-full.png)

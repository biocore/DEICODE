
[![Latest Github release](https://img.shields.io/github/release/cjm007/DEICODE.svg)](https://github.com/cjm007/DEICODE/releases/latest)
[![Build Status](https://travis-ci.org/cameronmartino/DEICODE.svg?branch=master)](https://travis-ci.org/cameronmartino/DEICODE)[![Coverage Status](https://coveralls.io/repos/github/cjm007/DEICODE/badge.svg?branch=master)](https://coveralls.io/github/cjm007/DEICODE?branch=master)
[![DOI](https://zenodo.org/badge/72654142.svg)](https://zenodo.org/badge/latestdoi/72654142)

deicode is a tool box for running Robust Aitchison RPCA on sparse omics datasets, linking specific features (i.e. microbes or metabolites) to beta-diversity ordination.

## Installation

To install the most up to date version of deicode, run the following command

    pip install git+https://github.com/cameronmartino/DEICODE.git

## Qiime2 tutorial

First make sure that qiime2 is installed before installing deicode. Then run

```
qiime dev refresh-cache
```

Once qiime2 is properly interfaced with deicode, you can import your biom tables
into Artifacts.  Here we will be using the [Sleep Apnea dataset](https://qiita.ucsd.edu/study/description/10422)
as an example. The full example run can be [found here](https://github.com/cameronmartino/DEICODE/blob/master/Examples/sleep_apnea/SleepApnea-qiime2-tutorial.ipynb)

```
qiime tools import \
    --input-path qiita_10422_table.biom \
    --output-path qiita_10422_table.biom.qza \
    --type FeatureTable[Frequency]
```
You can then run the qiime2 deicode rpca-biplot commmand as follows.

```
qiime deicode rpca-biplot \
    --i-table qiita_10422_table.biom.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot ordination.qza
```
Once you have this, you can directly visualize this in emperor
```
qiime emperor biplot \
    --i-biplot ordination.qza \
    --m-sample-metadata-file qiita_10422_metadata.tsv \
    --o-visualization biplot.qzv \
    --p-number-of-features 8
```
You can view the resulting visualization at https://view.qiime2.org.
It should look as follows
Inline-style:
![biplot](https://github.com/cameronmartino/DEICODE/blob/master/Examples/sleep_apnea/qiime_view.png)

## Usage

Python
```python

from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
import numpy as np

# rclr preprocessing

# numpy.ndarray - a array of counts (samples,features) with shape (M,N) where N>M
table=np.array([[3, 3, 0], [0, 4, 2], [3, 0, 1]]) 
table_rclr=rclr().fit_transform(table)

# OptSpace (RPCA)

opt=OptSpace().fit(table_rclr)
U=opt.sample_weights # numpy.ndarray - "Sample Loadings" 
V=opt.feature_weights # numpy.ndarray - "Feature Loadings" 
s=opt.s # numpy.ndarray - The singular values

# or directly transform

U,s,V=OptSpace().fit_transform(table_rclr)

```

Command line
```sh
Usage: deicode_rpca [OPTIONS]

  Runs RPCA with an rclr preprocessing step

Options:
  --in_biom TEXT              Input table in biom format.
  --output_dir TEXT           Location of output files.
  --rank INTEGER              Rank with witch to run OptSpace. default=3
  --min_sample_depth INTEGER  Minimum Sample Sequencing Depth Cut Off
                              default=500
  --help                      Show this message and exit.
```


## Simulation Benchmarking

* [simulations](https://github.com/cameronmartino/DEICODE/blob/master/Benchmarking/simulations/simulations.ipynb)
* [case studies](https://github.com/cameronmartino/DEICODE/blob/master/Benchmarking/case_studies/case_studies.ipynb)

## Other Resources

The code for OptSpace was translated to python from a [MATLAB package](http://swoh.web.engr.illinois.edu/software/optspace/code.html) maintained by Sewoong Oh (UIUC).

[Simulation and Case Study Benchmarking](https://github.com/cameronmartino/DEICODE/tree/master/Benchmarking)

- Transforms and PCoA : [Scikit-bio](https://github.com/biocore/scikit-bio)
- Data For Examples : [Qiita](https://qiita.ucsd.edu/)

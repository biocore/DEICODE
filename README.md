
[![Build Status](https://travis-ci.org/cameronmartino/DEICODE.svg?branch=master)](https://travis-ci.org/cameronmartino/DEICODE)
[![Coverage Status](https://coveralls.io/repos/github/cameronmartino/DEICODE/badge.svg)](https://coveralls.io/github/cameronmartino/DEICODE)

deicode is a tool box for running Robust Aitchison RPCA on sparse omics datasets, linking specific features to beta-diversity ordination.

Note that deicode is not compatible with python 2, and is compatible with Python 3.4 or later. deicode is currently in alpha. We are actively developing it, and backward-incompatible interface changes may arise.

## Installation

To install the most up to date version of deicode, run the following command

    # stable version (currently v.0.1.0)
    pip install deicode

    # dev. version
    pip install git+https://github.com/cameronmartino/DEICODE.git

## Tutorials

* [What is Robust Aitchison RPCA](https://github.com/cameronmartino/DEICODE/blob/master/ipynb/introduction.ipynb)

### Qiime2 tutorial

* [Sleep Apnea Biplots](https://github.com/cameronmartino/DEICODE/blob/master/ipynb/sleep_apnea/SleepApnea-qiime2-tutorial.ipynb)

First make sure that qiime2 is installed before installing deicode. Then run

```
qiime dev refresh-cache
```

Once qiime2 is properly interfaced with deicode, you can import your biom tables
into Artifacts.  Here we will be using the [Sleep Apnea dataset](https://qiita.ucsd.edu/study/description/10422)
as an example.

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
    --m-feature-metadata-file taxonomy.qza \
    --o-visualization biplot.qzv \
    --p-number-of-features 8
```
You can view the resulting visualization at https://view.qiime2.org.
It should look as follows
![biplot](https://github.com/cameronmartino/DEICODE/blob/master/ipynb/sleep_apnea/qiime_view.png)

### Python Tutorial

* [Sleep Apnea Log Ratio Tutorial](https://github.com/cameronmartino/DEICODE/blob/master/ipynb/sleep_apnea/SleepApnea-python-tutorial.ipynb)

## Simulation Benchmarking

* [simulations](https://github.com/cameronmartino/DEICODE/tree/master/benchmarking/simulations)
* [case studies](https://github.com/cameronmartino/DEICODE/tree/master/benchmarking/case_studies)

## Other Resources

The code for OptSpace was translated to python from a [MATLAB package](http://swoh.web.engr.illinois.edu/software/optspace/code.html) maintained by Sewoong Oh (UIUC).

[Simulation and Case Study Benchmarking](https://github.com/cameronmartino/DEICODE/tree/master/benchmarking)

- Transforms and PCoA : [Scikit-bio](https://github.com/biocore/scikit-bio)
- Data For Examples : [Qiita](https://qiita.ucsd.edu/)

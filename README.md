[![Build Status](https://travis-ci.org/biocore/DEICODE.svg?branch=master)](https://travis-ci.org/biocore/DEICODE)
[![Coverage Status](https://coveralls.io/repos/github/biocore/DEICODE/badge.svg?branch=master)](https://coveralls.io/github/biocore/DEICODE?branch=master)

DEICODE is a tool box for running Robust Aitchison PCA on sparse compositional omics datasets, linking specific features to beta-diversity ordination. 

## Installation

To install the most up to date version of deicode, run the following command

    # pip (only supported for QIIME2 >= 2018.8)
    pip install deicode

    # conda (only supported for QIIME2 >= 2019.1)
    conda install -c conda-forge deicode 

**Note**: that deicode is not compatible with python 2, and is compatible with Python 3.4 or later. deicode is currently in alpha. We are actively developing it, and backward-incompatible interface changes may arise.

## Using DEICODE as a standalone tool

```
$ deicode --help
Usage: deicode [OPTIONS]

  Runs RPCA with an rclr preprocessing step.

Options:
  --in-biom TEXT               Input table in biom format.  [required]
  --output-dir TEXT            Location of output files.  [required]
  --rank INTEGER               The underlying low-rank structure (suggested: 1
                               < rank < 10) [minimum 2]  [default: 3]
  --min-sample-count INTEGER   Minimum sum cutoff of sample across all
                               features  [default: 500]
  --min-feature-count INTEGER  Minimum sum cutoff of features across all
                               samples  [default: 10]
  --iterations INTEGER         The number of iterations to optimize the
                               solution (suggested to be below 100; beware of
                               overfitting) [minimum 1]  [default: 5]
  --help                       Show this message and exit.
```

## Using DEICODE inside [QIIME 2](https://qiime2.org/)

* The QIIME2 forum tutorial can be found [here](https://forum.qiime2.org/t/robust-aitchison-pca-beta-diversity-with-deicode/8333).
* The official plugin docs and tutorial can be found [here](https://library.qiime2.org/plugins/deicode).
* The in-repo tutorial can be found [here](https://github.com/biocore/DEICODE/blob/master/ipynb/tutorials/moving-pictures.md).

## Other Resources

* [Aitchison Distance Introduction](https://github.com/biocore/DEICODE/blob/master/ipynb/introduction.ipynb)

- The code for OptSpace was translated to python from a [MATLAB package](http://swoh.web.engr.illinois.edu/software/optspace/code.html) maintained by Sewoong Oh (UIUC).
- Transforms and PCoA : [Scikit-bio](http://scikit-bio.org)
- Data For Examples : [Qiita](https://qiita.ucsd.edu/)

#### Simulation and Benchmarking

* [simulations and case studies](https://github.com/cameronmartino/deicode-benchmarking)

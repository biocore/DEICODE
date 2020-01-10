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

## Using DEICODE inside [QIIME 2](https://qiime2.org/)

* The QIIME2 forum tutorial can be found [here](https://forum.qiime2.org/t/robust-aitchison-pca-beta-diversity-with-deicode/8333).
* The official plugin docs and tutorial can be found [here](https://library.qiime2.org/plugins/deicode).
* The in-repo tutorial can be found [here](https://nbviewer.jupyter.org/github/biocore/DEICODE/blob/master/ipynb/tutorials/moving-pictures.ipynb).

## Using DEICODE as a standalone tool

There are two commands within deicode. The first is `rpca` and the second is `auto-rpca`. The only difference is that `auto-rpca` automatically estimates the underlying-rank of the matrix and requires no input for the `n_components` parameter. In the `rpca` the `n_components` must be set explicitly. The structure of the commands follows the QIIME2 commands exactly and so questions about the use of the tool can be answered in the tutorial in the `Using DEICODE inside QIIME 2` section above. However, an example analysis without the use of QIIME2 can be found [here](https://nbviewer.jupyter.org/github/biocore/DEICODE/blob/master/ipynb/tutorials/moving-pictures-standalone-cli-and-api.ipynb).

## Using DEICODE as a Python API

The `rpca` functionality of DEICODE is also exposed as a python API. An example analysis without the use of the command line can be found [here](https://nbviewer.jupyter.org/github/biocore/DEICODE/blob/master/ipynb/tutorials/moving-pictures-standalone-cli-and-api.ipynb).

## Other Resources

* [Aitchison Distance Introduction](https://github.com/biocore/DEICODE/blob/master/ipynb/introduction.ipynb)

- The code for OptSpace was translated to python from a [MATLAB package](http://swoh.web.engr.illinois.edu/software/optspace/code.html) maintained by Sewoong Oh (UIUC).
- Transforms and PCoA : [Scikit-bio](http://scikit-bio.org)
- Data For Examples : [Qiita](https://qiita.ucsd.edu/)

#### Simulation and Benchmarking

* [simulations and case studies](https://github.com/cameronmartino/deicode-benchmarking)

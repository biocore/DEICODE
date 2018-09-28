
[![Latest Github release](https://img.shields.io/github/release/cjm007/DEICODE.svg)](https://github.com/cjm007/DEICODE/releases/latest)
[![Build Status](https://travis-ci.org/cameronmartino/DEICODE.svg?branch=master)](https://travis-ci.org/cameronmartino/DEICODE)[![Coverage Status](https://coveralls.io/repos/github/cjm007/DEICODE/badge.svg?branch=master)](https://coveralls.io/github/cjm007/DEICODE?branch=master)
[![DOI](https://zenodo.org/badge/72654142.svg)](https://zenodo.org/badge/latestdoi/72654142)

deicode is a tool box for running Robust Aitchison RPCA on sparse omics datasets, linking specific features (i.e. microbes or metabolites) to beta-diversity ordination.

## Installation
    
To install the most up to date version of DEICODE, run the following command

    pip install git+https://github.com/cameronmartino/DEICODE.git

## Examples

[Sponge Life Stage Case Study](https://github.com/cameronmartino/DEICODE/blob/master/Examples/sponge_biom.ipynb)

[Infant Development Study](https://github.com/cameronmartino/DEICODE/blob/master/Examples/infant_gut.ipynb)

## Usage

    from deicode.optspace import OptSpace
    from deicode.preprocessing import rclr
    import numpy as np

    # rclr preprocessing

    table_norm=rclr().fit_transform(otutabledf.copy())

    # OptSpace (RPCA)

    opt=OptSpace().fit(table_norm)
    U=opt.sample_weights # numpy.ndarray - "Sample Loadings" 
    V=opt.feature_weights # numpy.ndarray - "Feature Loadings" 
    s=opt.s # numpy.ndarray - The singular values
    result=opt.solution # numpy.ndarray - (U*S*V.transpose()) of shape (M,N)
    
    # or 

    U,s,V=OptSpace().fit_transform(table_norm)
    result=np.dot(np.dot(U,s),v.T)

## Simulation Benchmarking 

[rclr preprocessing](https://github.com/cameronmartino/DEICODE/blob/master/Benchmarking/simulations.ipynb)

## Other Resources 

The code for OptSpace was translated to python [MATLAB package](http://swoh.web.engr.illinois.edu/software/optspace/code.html) maintained by Sewoong Oh (UIUC).

- Transforms and PCoA : [Scikit-bio](https://github.com/biocore/scikit-bio)
- Other Imputation Methods: [Fancy Impute](https://github.com/hammerlab/fancyimpute)
- Data For Examples : [Qiita](https://qiita.ucsd.edu/)

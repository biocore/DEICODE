[![Build Status](https: // travis - ci.org / biocore / DEICODE.svg?branch=master)](https: // travis - ci.org / biocore / DEICODE)
[![Coverage Status](https: // coveralls.io / repos / github / biocore / DEICODE / badge.svg?branch=master)](https: // coveralls.io / github / biocore / DEICODE?branch=master)

Deicode is a tool box for running Robust Aitchison PCA on sparse compositional omics datasets, linking specific features to beta - diversity ordination.

# Installation

To install the most up to date version of deicode, run the following command

# pip (only supported for Qiime >= 2018.8)
pip install deicode

# conda (only supported for Qiime >= 2019.1)
conda install - c conda - forge deicode

**Note**: that deicode is not compatible with python 2, and is compatible with Python 3.4 or later. deicode is currently in alpha. We are actively developing it, and backward - incompatible interface changes may arise.

# Qiime2 tutorial

**Note**: the official plugin docs and tutorial can be found[here](https: // library.qiime2.org / plugins / q2 - deicode).

* [Sleep Apnea Biplots](https: // github.com / biocore / DEICODE / blob / master / ipynb / sleep_apnea / SleepApnea - qiime2 - tutorial.md)

# Other Resources

* [What is Robust Aitchison RPCA](https: // github.com / biocore / DEICODE / blob / master / ipynb / introduction.ipynb)

- The code for OptSpace was translated to python from a[MATLAB package](http: // swoh.web.engr.illinois.edu / software / optspace / code.html) maintained by Sewoong Oh(UIUC).
- Transforms and PCoA: [Scikit - bio](http: // scikit - bio.org)
- Data For Examples: [Qiita](https: // qiita.ucsd.edu /)

#### Simulation and Benchmarking

* [simulations](https: // github.com / biocore / DEICODE / tree / master / benchmarking / simulations)
* [case studies](https: // github.com / biocore / DEICODE / tree / master / benchmarking / case_studies)

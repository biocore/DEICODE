# DEICODE
## Discovery of Environmental Influences through Convex Optimized Decomposition by Ecotypes (DEICODE) 

Through matrix completion methods trends can be discovered from environmental changes between sample groups.

This project is a demonstration of this method in 16S rRNA sequencing data. 

## Matrix Completion Examples

### [88 Soils DEICODE Method Example](https://github.com/cjm007/DEICODE/blob/master/88Soils.ipynb)
### [Crohns DEICODE Method Example](https://github.com/cjm007/DEICODE/blob/master/Crohn_example.ipynb)

## Benchmarking Compltion Methods

### [High-rank Challange 88 Soils Imputation Benchmarking](https://github.com/cjm007/DEICODE/blob/master/Imputation_88Soils_benchmarking.ipynb)

## Installation

cd usr../DEICODE

chmod +x setup.sh

./setup.sh

## Test

cd usr../DEICODE

source activate DEICODE_env

./DEICODE.py -i data/88_soils.biom  -m data/soil_map.txt -o test_output

## Usage

cd usr.../DEICODE 

source activate DEICODE_env 

    usage: DEICODE.py [-h] -i INPUT_OTU -m MAP -o OUTPUT [-l LOW_RANK_METHOD]
                      [-d DECOMPIT] [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE]
                      [-s MAPSTART] [-f FEATURE] [-w_zero W_ZERO] [-w_low W_LOW]
                      [-w_high W_HIGH] [-n NCOMP] [-fm FMETHOD]


### Commands 

    usage: DEICODE.py [-h] -i INPUT_OTU -m MAP -o OUTPUT [-l LOW_RANK_METHOD]
                      [-d DECOMPIT] [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE]
                      [-s MAPSTART] [-f FEATURE] [-w_zero W_ZERO] [-w_low W_LOW]
                      [-w_high W_HIGH] [-n NCOMP] [-fm FMETHOD]


Required arguments:

      -i INPUT_OTU, --Input_OTU INPUT_OTU
                            Path to .biom table i.e. home/usr/input/otu.biom (taxa
                            should be included in metadata)

      -m MAP, --map MAP     Path to Qiime style metadata i.e.
                            home/usr/input/map.txt

      -o OUTPUT, --output OUTPUT
                            Output directory

Optional arguments:

    -l LOW_RANK_METHOD, --low_rank_method LOW_RANK_METHOD
                          Specify a low rank method to use (default is
                          SoftImpute) (options = NNM (Nuclear Norm
                          Minimization), SoftImpute,IterativeSVD,
                          MatrixFactorization, can also use WPCA, or EMPCA for
                          imputation)

    -d DECOMPIT, --decompit DECOMPIT
                          How many iterations to complete in decomposition
                          (default=100) (options = any integer or None)

    -b BACTNUM, --bactnum BACTNUM
                          Number of bacteria to extract from PCA axis
                          (default=12) (options = any integer)

    -c CLASSNUM, --classnum CLASSNUM
                          Number of highest scoring classifiers to use in
                          analysis (default=2) (options = any integer greater
                          than 1 and less than the umber of columns in the
                          mapping file)

    -t TAXAUSE, --taxause TAXAUSE
                          What level of taxonomy to extract from PCA axis
                          (default=genus) (options = phylum, class, order,
                          family, genus, species or None if you do not have
                          incorporated taxaonomy)

    -s MAPSTART, --mapstart MAPSTART
                          What column to start analysis on in mapping file,
                          (i.e. skipping barcode sequences) (default=3) (options
                          = any integer greater than 1 and less than the umber
                          of columns in the mapping file)

    -f FEATURE, --feature FEATURE
                          Set to False if you would like to turn off feature
                          selection (default=True, warning: on large datastes
                          this will signficantly slow down run time)

    -w_zero W_ZERO, --w_zero W_ZERO
                          Zero value used for weighted PCA (default=0)

    -w_low W_LOW, --w_low W_LOW
                          Low value used for weighted PCA (default .001)

    -w_high W_HIGH, --w_high W_HIGH
                          High value used for weighted PCA (default 10)

    -n NCOMP, --ncomp NCOMP
                          Number of Principle Components to Use in Feature
                          Selection (default=3)

    -fm FMETHOD, --fmethod FMETHOD
                          Method to use for feature selection (default=WPCA)

## Credits


- Source Imputation: https://github.com/hammerlab/fancyimpute
- Source WPCA/EMPCA: https://github.com/jakevdp/wpca
- Source decompostion Soil and Crohn example: https://github.com/dganguli/robust-pca

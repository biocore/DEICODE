# DEICODE
## Discovery of Environmental Influences through Convex Optimized Decomposition by Ecotypes (DEICODE) 

Through matrix completion methods trends can be discovered from environmental changes between sample groups.

This project is a demonstration of this method in 16S rRNA sequencing data. 

## Examples Using Nuclear Norm Rank Minimization

### 88 Soils - https://github.com/cjm007/DEICODE/blob/master/88Soils.ipynb
### Crohns - https://github.com/cjm007/DEICODE/blob/master/Crohn_example.ipynb

## Benchmarking Examples

### 88 Soils - coming soon
### Crohns - coming soon

## Installation

cd usr../DEICODE

chmod +x setup.sh

./setup.sh

## Usage

cd usr.../DEICODE 

source activate DEICODE_env 

./DEICODE.py [-h] -i INPUT_OTU -m MAP -o OUTPUT [-l LOW_RANK_METHOD]
                  [-d DECOMPIT] [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE]
                  [-s MAPSTART]

### Commands 

usage: DEICODE.py [-h] -i INPUT_OTU -m MAP -o OUTPUT [-l LOW_RANK_METHOD]
                  [-d DECOMPIT] [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE]
                  [-s MAPSTART]

 required arguments:
 
  -i INPUT_OTU, --Input_OTU INPUT_OTU
                        Path to .biom table i.e. home/usr/input/otu.biom (taxa
                        should be included in metadata)
                        
  -m MAP, --map MAP     Path to Qiime style metadata i.e.
                        home/usr/input/map.txt
                        
  -o OUTPUT, --output OUTPUT

                        
 optional arguments:
 
  -l LOW_RANK_METHOD, --low_rank_method LOW_RANK_METHOD
                        Specify a low rank method to use (default is
                        IterativeSVD) (options = NNM (Nuclear Norm
                        Minimization), SoftImpute,IterativeSVD, MICE, KNN,
                        WPCA, or EMPCA)
                        
  -d DECOMPIT, --decompit DECOMPIT
                        How many iterations to complete in decomposition
                        (deafault=100) (options = any integer)
                        
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
                        (deafult=genus) (options = phylum, class, order,
                        family, genus, species or None if you do not have
                        incorporated taxaonomy)
                        
  -s MAPSTART, --mapstart MAPSTART
                        What column to start analysis on in mapping file,
                        (i.e. skipping barcode sequences) (deafult=3) (options
                        = any integer greater than 1 and less than the umber
                        of columns in the mapping file)
## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

Source decomposition code and algorithm is from https://github.com/dganguli/robust-pca

Soil and Human microbiome examples of use in biological data written by cameron martino 

### email cameronmartino at gmail dot com for concerns or questions


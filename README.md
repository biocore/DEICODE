# DEICODE
## Discovery of Environmental Influences through Convex Optimized Decomposition by Ecotypes (DEIECODE) 


By decomposing the sparse data into its dense low-rank component trends can be discorved from enviornmental cgnages between sample groups.

This project is a demonstration of this method in 16S rRNA sequencing data. 

## Installation

cd usr../DEIECODE

chmod +x setup.sh

./setup.sh

## Examples

cd usr.../DEIECODE 

source activate COD_env 

jupyter notebook 

A browser window will open, then open notebooks to view examples from paper. Shift-enter to run cells. 


## Usage

cd usr.../DEIECODE 

source activate COD_env 

./DEICODE.py [-h] -i INPUT_DIR -m MAP -o OUTPUT [-d DECOMPIT]
                   [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE] [-s MAPSTART]


### Commands 

usage: DEIECODE.py [-h] -i INPUT_DIR -m MAP -o OUTPUT [-d DECOMPIT]
                   [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE] [-s MAPSTART]

required arguments:

  -i INPUT_DIR, --input_dir INPUT_DIR
                        path to .biom table i.e. home/usr/input/otu.biom

  -m MAP, --map MAP     path to Qiime style metadata i.e.
                        home/usr/input/map.txt

  -o OUTPUT, --output OUTPUT
                        Output directory

optional arguments:

  -h, --help            show this help message and exit

  -d DECOMPIT, --decompit DECOMPIT
                        How many iterations to complete in decomposition
                        (deafault=48) (options = any integer)

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
                        (deafult=family) (options = phylum, class, order,
                        family, genus, species)

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

for issues please contact cameronmartino at gmail dot com 


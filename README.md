# DEICODE
## Discovery of Environmental Influences through Convex Optimized Decomposition by Ecotypes (DEICODE) 

By decomposing the sparse data into its dense low-rank component trends can be discorved from enviornmental cgnages between sample groups.

This project is a demonstration of this method in 16S rRNA sequencing data. 

## Introduction to matrix decomposition and motivation for its application to 16S analysis 

#### What is PCA
* Standardize the data.
* Obtain the Eigenvectors and Eigenvalues from the covariance matrix or correlation matrix, or perform Singular Vector Decomposition.
* Sort eigenvalues in descending order and choose the k eigenvectors that correspond to the k largest eigenvalues where k is the number of dimensions of the new feature subspace (k < d).
* Construct the projection matrix W from the selected k eigenvectors.
* Transform the original dataset X via W to obtain a k-dimensional feature subspace Y.

The eigenvectors and eigenvalues of a covariance (or correlation) matrix represent the "core" of a PCA: The eigenvectors (principal components) determine the directions of the new feature space, and the eigenvalues determine their magnitude. In other words, the eigenvalues explain the variance of the data along the new feature axes.

#### Why would you use PCA in 16S analysis 

* Dimentionality Reduction and visualization in 3D

16S rRNA analysis begins by processing raw seaquencing reads into a matrix called and OTU table that is high dimentional and by using PCA the data can be visaluzed in 3 dimentions or less.

* Extract what bacteria are cauing variance between samples directly from PCA visualization 

From PCA graphs you can view which bacteria are contributing most to the chnages along an axis in your graph. This means unlike the conventional PCoA, in PCA your axis has meaning!

* Support vector machines (SVMs) can be used to quickly determine the best column from metadata (i.e pH from map.txt)

If you can perfrom PCA is can be used as input for SVMs, a very poerful tool as metadata continues to grow in size. 

#### Why has no one used this before?

##### There are lots of reasons (more than what is listed)

* OTU tables consist of m rows of OTUs each representing a potential microbe and n columns of samples where m is much greater than n. 
* OTU tables are very sparse, meaning they contain alot of zeros 
* Compositionality sum constraints skew variance-covariance measurements preventing the use of multivariate analysis that rely on multivariate normality. 

#### The Solution Convex Optomized Decomposition!

By using the l1 and nuclear norm we can decompose the matrix into its low-rank and sparse compoenents and then use the low-rank matrix for PCA. This is simillar to reucing noise in an image, allowing you to see the trends.

i.e. 

![alt tag](https://github.com/cjm007/DEICODE/blob/master/etc/decomp.png)

Where we can used the picture in the middle to determine trends, here the rank is three and the pattern is a checker board with a grey and white side. 

This allows us to visualize the data using PCA, use SVM to determine the best classifier from metadata, and we can extract what bacteria are cuasing variance between your groups determined by SVM. 

### Check out the example jupyter notebooks above by just opening them in your browser on git hub.  email cameronmartino at gmail dot com for concerns or questions

#
#
## Installation

cd usr../DEICODE

chmod +x setup.sh

./setup.sh

## Examples

cd usr.../DEICODE 

source activate DEICODE_env 

jupyter notebook 

A browser window will open, then open notebooks to view examples from paper. Shift-enter to run cells. 


## Usage

cd usr.../DEICODE 

source activate DEICODE_env 

./DEICODE.py [-h] -i INPUT_DIR -m MAP -o OUTPUT [-d DECOMPIT]
                   [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE] [-s MAPSTART]


### Commands 

usage: DEIECODE.py [-h] -i INPUT_DIR -m MAP -o OUTPUT [-d DECOMPIT]
                   [-b BACTNUM] [-c CLASSNUM] [-t TAXAUSE] [-s MAPSTART]

-i Path to a .biom file with taxanomic metadata 

-m Path to Qiime style tab delimeted metadata (http://qiime.org/documentation/file_formats.html)

-o Output directory for visualization and csv file results 


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


#!/bin/bash

# cameronmartino@gmail.com for bugs

dir=$(pwd)
conda env create -f environment.yml
git clone https://github.com/dganguli/robust-pca.git $dir/r_pca
chmod +x DEICODE.py

echo Done! Please activate the envrionment DEICODE_env and navigate to the jupyter notebooks
echo To source activate the envrionment execute:$ source activate DEICODE_env

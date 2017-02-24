#!/usr/bin/env bash

# cameronmartino@gmail.com for bugs


dir=$(pwd)

if [ "$(uname)" == "Linux" ]; then
conda env create -f environment_linux.yml
elif [ "$(uname)" == "Darwin" ]; then
conda env create -f environment.yml
else
echo Will only run on linux and OSX
fi

chmod +x  DEICODE.py


echo Done! Please activate the envrionment DEICODE_env and navigate to the jupyter notebooks
echo To source activate the envrionment execute:$ source activate DEICODE_env

#!/bin/bash

# untar your Python installation
tar -xzf parallel_grid_search.tar.gz
# make sure the script will use your Python installation,
# and the working directory as it's home location
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home

# buy some time
# sleep 5m

# run your script, $1 will be the line in the file (i.e., params dict)
python param_optimiz_RFC_multilabel_full.py $1

# clean up
rm 50percentplusTraining.csv param_optimiz_RFC_multilabel.py RFC_params.tsv

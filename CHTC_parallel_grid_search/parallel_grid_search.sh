#!/bin/bash

# untar your Python installation
tar -xzf parallel_grid_search.tar.gz
# make sure the script will use your Python installation,
# and the working directory as it's home location
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home

# buy some time
sleep 5m

# run your script, $1 will be the line in the file (i.e., params dict)
python cross_val_score.py $1

# clean up
rm Features_AvgMaxCV_Weights.csv cross_val_score.py params.tsv

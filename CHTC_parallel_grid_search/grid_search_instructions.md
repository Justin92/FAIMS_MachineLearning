#Grid Search Instructions
Specific file labels pertain to the test grid search using multilabel random forest classifier,
but can easily be made widely applicable
##1 Write cross validation Code

Currently known functional cross validation code in
Dekstop/FAIMS_MachineLearning/CHTC_parallel_grid_search/param_optimiz_RFC_multilabel.py

Some things to keep in mind about cross validation code:
* If there are floats in the parameter table, it will often interperet ints as floats so
wrap variables that need to be integers in int()
* Try to do as much of preprocessing  outside of code as possible and get rid of
excess code snippets CHTC python will run into errors if code indented below comment
(e.g body of for loop not fully commmented out)

##2 Generate parameters index and table

Currently two python modules for making parameter tables in
Dekstop/FAIMS_MachineLearning/CHTC_parallel_grid_search/

**make_params_list.py** and **make_params_short_list.py** modify the latter to
generate the params you want

```bash
  rm short_params.tsv  short_params_indices                  #Remove old param files
  python make_params_short_list.py                 #Generate new params files
  nano short_params_indices                        #Go into indices and delete
```

Then you will need to go into and delete th column header "index" in the short_params_indices


##3 Create correct python environment and files on to CHTC (4 files)

Will then need to use scp to get your data csv, short_params_indices, short_params.tsv
and param_optimiz_RFC_multilabel.py (your cross validation code) onto CHTC

If you're in the correct directory:

```bash
  scp /p/JGM_FAIMS_CVprediction/JMM_PreProcessed_Data/50percentplusTraining.csv mcketney@submit-1.chtc.wisc.edu:/home/mcketney/CHTC_multilabel_grid_search/
  scp short_params_indices mcketney@submit-1.chtc.wisc.edu:/home/mcketney/CHTC_multilabel_grid_search/
  scp short_params.tsv mcketney@submit-1.chtc.wisc.edu:/home/mcketney/CHTC_multilabel_grid_search/
  scp param_optimiz_RFC_multilabel.py mcketney@submit-1.chtc.wisc.edu:/home/mcketney/CHTC_multilabel_grid_search/

```

Once files moved will need to make sure correct references in parallel_grid_search.sub,
parallel_grid_search.sh and param_optimiz_RFC_multilabel.py (line 65-80)

Will also need to build correct python installation using interactive node

If starting on local computer(based on [Ian's repo](https://github.com/ijmiller2/CHTC-RFR/blob/master/README.md))

```bash
  ssh mcketney@submit-1.chtc.wisc.edu
    #enter password
  cd CHTC_multilabel_grid_search/
  condor_submit -i interactive_compile.sub  #Comes from Ians repository
    #greeting welcoming you to node
  wget https://www.python.org/ftp/python/3.6.8/Python-3.6.8.tgz
  tar zxvf Python-3.6.8.tgz
  mkdir python
  cd Python-3.6.8
  ./configure --prefix=$(pwd)/../python
  make
  make install
  cd ..
  cp python/bin/python3 python/bin/python
  cp python/bin/pip3 python/bin/pip
  export PATH=$(pwd)/python/bin:$PATH
  which pip
  pip install http://download.pytorch.org/whl/cpu/torch-1.0.0-cp36-cp36m-linux_x86_64.whl
  pip install fastai
  pip install scikit-multilearn
  pip install  -U imbalanced-learn
  pip install pandas
  pip install sklearn
  pip install xgboost
  tar -czvf python.tar.gz python/
  exit
  #Now you're back in your submit node
```

```bash
  tar zxvf python.tar.gz   #Untar your python install into python/
  #Repackage python, data and param_table into tar file
  tar zcvf parallel_grid_search.tar.gz python param_optimiz_RFC_multilabel.py 50percentplusTraining.csv short_params.tsv

  #submit job to condor
  condor_submit parallel_grid_search.sub
```

##4 Concatenate files and get down from CHTC

Find the easiest strategy is to move all tsv files to single directory then just
move back the ones I don't want

```bash
mkdir tsv_outputs/
mv *.tsv tsv_outputs/
cp concatenate_output_dfs.py
mv tsv_outputs/short_params.tsv /home/mcketney/CHTC_multilabel_grid_search/
cd tsv_ouptuts/
python concatenate_output_dfs.py
```
To move to pull files down from CHTC you have to exit to your local computer

```bash
  exit
  scp mcketney@submit-1.chtc.wisc.edu:/home/mcketney/CHTC_multilabel_grid_search/tsv_outputs/model_performance.tsv /c/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/CHTC_parallel_grid_search/multilabel_grid_search/
```
That's that

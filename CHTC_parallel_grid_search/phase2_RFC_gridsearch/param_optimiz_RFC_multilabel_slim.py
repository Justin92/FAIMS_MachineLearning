
import sys
import timeit
import numpy as np
import pandas as pd
import sklearn
import sklearn.ensemble
#import keras
import fastai
import skmultilearn

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, GridSearchCV
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.dummy import DummyRegressor
from sklearn import ensemble
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score
from numpy import array
from numpy import argmax                   #finds the index of the maximum value in a vector
#from keras.utils import to_categorical
#from keras_preprocessing.sequence import pad_sequences
from sklearn.multiclass  import OneVsRestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from imblearn.metrics import geometric_mean_score
#Gonna need to install the skmultilearn to do this
from skmultilearn.model_selection import IterativeStratification
#Gonna need to install fastai for this
from fastai.metrics import fbeta
from fastai.torch_core import np2model_tensor
from fastai.metrics import accuracy_thresh




#####READING IN DATA, ESTABLISHING X AND Y#########
#Read in data that already has all the pyteomics features added
#and potentially one-hot
data_df = pd.read_csv("50percentplusTraining.csv", low_memory=False) #P:\JGM_FAIMS_CVprediction\JMM_PreProcessed_Data\50percentplusTraining.csv

#Ideally input data structured with labels, followed by named features, followed
#by one-hots
onehotstart = 22

#Make X from concatenating named features to one-hots
feature_subset = ['Charge', 'Length', 'pyMass', 'pI']
X = np.concatenate((data_df[feature_subset], data_df.iloc[:, onehotstart:]), axis = 1)

# Really want to try more liberal threshold for training but think that will be
# difficult here
y = data_df.loc[ : ,  'X20':'X95'].values




#########################RANDOM FOREST REGRESSOR START##################################



####INSTANTIATING AND SCORING THE MODEL########


### run single parameter set (using stdin) ###
params_index = int(sys.argv[1])
params_df = pd.read_csv("RFC_params2.tsv", sep = "\t")
# grab the parameter set based on the index
params = params_df.iloc[params_index]
# initial the dataframe to record model parameters and performance
score_df = pd.DataFrame([params], index=[params_index])
multiforest = RandomForestClassifier(   #May be easier to control OVR RFC

    max_depth=params['max_depth'],
    #min_samples_leaf=int(params['min_samples_leaf']),
    min_samples_split=int(params['min_samples_split']),
    n_estimators=int(params['n_estimators']),
    max_features=params['max_features'],  #This is a new one.
    #Don't remember if strings and numbers can work
    random_state=23
)
#)
# get the cross validation scores, keep track of time
start = timeit.default_timer()

RFC_geomean = []
RFC_f2 = []
RFC_precision = []
RFC_recall = []
RFC_rocauc = []
RFC_threshacc = []

k_fold = IterativeStratification(n_splits=5, order=3)#, random_state=123)
for train, test in k_fold.split(X, y):

        #Train the multiforest using the training indices
        multiforest.fit(X[train], y[train])

        #Generate predictions
        y_predict = multiforest.predict(X[test])
        #Generate probs if want to use RFC OVR
        #y_predprobs = multiforest.predict_proba(X[test])
        #If want to use innate multilabel capacity of RFC
        y_problist = multiforest.predict_proba(X[test])
        y_predprobs = np.array([[]] *X[test].shape[0])
        for item in y_problist:
            y_predprobs = np.concatenate((y_predprobs, np.reshape(item[:, 1], (-1, 1))), axis = 1)

        #Gotta convert both of them to tensors hahahaha
        RFC_f2.append(fbeta(np2model_tensor(y_predprobs), np2model_tensor(y[test]), thresh=0.55, beta =2, sigmoid=False))
        RFC_threshacc. append(accuracy_thresh(np2model_tensor(y_predprobs), np2model_tensor(y[test]), thresh=0.55, sigmoid=False))
        RFC_precision.append(average_precision_score( y[test], y_predprobs, average='weighted'))
        RFC_rocauc.append(roc_auc_score( y[test], y_predprobs, average='weighted'))
        RFC_recall.append(recall_score( y[test], y_predict, average='weighted'))
        gmean = 0
        for column in range(0,y_predict.shape[1]):
            gmean = gmean + geometric_mean_score(y[test, column], y_predict[:, column])
        RFC_geomean.append(gmean/y_predict.shape[1])


stop = timeit.default_timer()
# add the mean
score_df['Avg_Precision'] = round(np.mean(RFC_precision),6)
score_df['Avg_Recall'] = round(np.mean(RFC_recall),6)
score_df['Avg_ROCAUC'] = round(np.mean(RFC_rocauc),6)
score_df['Avg_ThrshAcc'] = round(np.mean(RFC_threshacc),6)
score_df['Avg_F2'] = round(np.mean(RFC_f2),6)
score_df['Avg_GeoMean'] = round(np.mean(RFC_geomean),3)




# get the runtime
run_time = stop - start
score_df['run_time'] = round(run_time,3)

# write out the file
outfile_name = "{}.tsv".format(params_index)
score_df.to_csv(outfile_name, index=False, sep = "\t") # index already built into table

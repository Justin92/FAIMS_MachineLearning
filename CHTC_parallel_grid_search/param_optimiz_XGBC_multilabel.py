
import sys
import timeit
import numpy as np
import pandas as pd
import sklearn
import sklearn.ensemble
import keras
import xgboost

from sklearn.ensemble import RandomForestRegressor
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
from keras.utils import to_categorical
from keras_preprocessing.sequence import pad_sequences
from sklearn.multiclass  import OneVsRestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from imblearn.metrics import geometric_mean_score
#Gonna need to install the skmultilearn to do this
from skmultilearn.model_selection import IterativeStratification
#Gonna need to install fastai for this
from fastai import fbeta
from fastai import np2model_tensor
from fastai import accuracy_thresh
# Gonna need to install XGBoost
from xgboost import XGBClassifier



#####READING IN DATA, ESTABLISHING X AND Y#########
#Read in data that already has all the pyteomics features added
#and potentially one-hot
data_df = pd.read_csv("FinalThreshold.csv", low_memory=False)

#Ideally input data structured with labels, followed by named features, followed
#by one-hots
onehotstart = #40  not sure what this is right now

#Make X from concatenating named features to one-hots
feature_subset = ['Charge', 'Length', 'pyMass', 'pI']
X = np.concatenate((data_df[feature_subset], data_df.iloc[:, onehotstart:], axis = 1)

# Really want to try more liberal threshold for training but think that will be
# difficult here
y = data_df.loc[ : ,  'X20':'X95'].values




#########################RANDOM FOREST REGRESSOR START##################################





###This section could be useful but I think it's faster
###to have the one-hot attached beforehand
#Adding one hot
#Start by finding the max and calculating needed vector length
VEC_LENGTH = max(evidence_df['Length']) * 20
#Define what residues are possible
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY '
#TURNING CHARACTERS INTO INTEGERS
# Map character keys to integer values in a dictionary, then map integer keys to character values to revers transform
char_to_int = dict((c, i) for i, c in enumerate(AMINO_ACIDS))   #character keys to integer values
int_to_char = dict((i, c) for i, c in enumerate(AMINO_ACIDS))   #integer keys to character values
#will store the sequence in a list
hotlist = list()
#Iterate through all rows of evidence
for i in range(0, evidence_df.shape[0]):
    pep = evidence_df['Sequence'][i]
    integer_encode = [char_to_int[char] for char in pep]
    encoded = to_categorical(integer_encode, num_classes=20)
    flatencode = encoded.flatten()
    #Add the flatened one hot sequence to my list
    hotlist.append(flatencode)
#Pad the ragged end of all the items in list so that have same number of columns
padded = pad_sequences(hotlist, padding= 'post', maxlen=VEC_LENGTH)
#Make that list into an array
hotarray = np.array(padded)


####INSTANTIATING AND SCORING THE MODEL########

RFC_geomean = []
RFC_f2 = []
RFC_precision = []
RFC_recall = []
RFC_rocauc = []
RFC_threshacc = []

### run single parameter set (using stdin) ###
params_index = int(sys.argv[1])
params_df = pd.read_csv("params.tsv", sep = "\t")
# grab the parameter set based on the index
params = params_df.iloc[params_index]
# initial the dataframe to record model parameters and performance
score_df = pd.DataFrame([params], index=[params_index])
xgbclf = OneVsRestClassifier(        #May be easier to control OVR RFC
    XGBClassifier(
    max_depth=params['max_depth'],
    learning_rate=params['learning_rate'],
    booster=params['booster'],
    max_delta_step=params['max_delta_step'],
    scale_pos_weight=params['scale_pos_weight']  #This is a new one.
    #Don't remember if strings and numbers can work
    random_state=23
)
)
# get the cross validation scores, keep track of time
start = timeit.default_timer()

XGB_geomean = []
XGB_f2 = []
XGB_precision = []
XGB_recall = []
XGB_rocauc = []
XGB_threshacc = []

k_fold = IterativeStratification(n_splits=5, order=4)#, random_state=123)
for train, test in k_fold.split(X, y):

        #Train the multiforest using the training indices
        xgbclf.fit(X[train], y[train])

        #Generate predictions
        y_predict = xgbclf.predict(X[test])
        #Generate probs if want to use RFC OVR
        y_predprobs = xgbclf.predict_proba(X[test])


        #Gotta convert both of them to tensors hahahaha
        XGB_f2.append(accuracy_thresh(np2model_tensor(y_predprobs), np2model_tensor(y1[test]), thresh=0.55, beta =2, sigmoid=False))
        XGB_threshacc. append(accuracy_thresh(np2model_tensor(y_predprobs), np2model_tensor(y1[test]), thresh=0.55, sigmoid=False))
        XGB_precision.append(average_precision_score( y1[test], y_predict, average='weighted'))
        XGB_rocauc.append(roc_auc_score( y1[test], y_predprobs, average='weighted'))
        XGB_recall.append(recall_score( y1[test], y_predict, average='weighted'))
        gmean = 0
        for(column in range(0,y_predict.shape[1])):
            gmean = gmean + geometric_mean_score(y1[test, column], y_predict[:, column])
        XGB_geomean.append(gmean/y_predict.shape[1])


stop = timeit.default_timer()
# add the mean
score_df['Avg_Precision'] = round(np.mean(XGB_precision),6)
score_df['Avg_Recall'] = round(np.mean(XGB_recall),6)
score_df['Avg_ROCAUC'] = round(np.mean(XGB_rocauc),6)
score_df['Avg_ThrshAcc'] = round(np.mean(XGB_threshacc),6)
score_df['Avg_F2'] = round(np.mean(XGB_f2),6)
score_df['Avg_GeoMean'] = round(np.mean(XGB_geomean),3)




# get the runtime
run_time = stop - start
score_df['run_time'] = round(run_time,3)

# write out the file
outfile_name = "{}.tsv".format(params_index)
score_df.to_csv(outfile_name, index=False, sep = "\t") # index already built into table

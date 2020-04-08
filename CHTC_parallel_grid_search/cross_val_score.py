
import sys
import timeit
import numpy as np
import pandas as pd
import sklearn
import sklearn.ensemble
import keras

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


#Read in data that already has all the pyteomics features added
evidence_df = pd.read_csv("Features_AvgMaxCV_Weights.csv", low_memory=False)

####FILTERING######
####Filter based on PEP
  #makes vector of indices of all PEPs above
badscores = evidence_df[evidence_df['PEP'] > 0.01].index
#Delete these bad scoring indices from the array
evidence_df.drop(badscores, inplace=True)
#Need to reset the indices of your table so they match much like rownames in R
evidence_df.index = range(len(evidence_df))

# downsample evidence_df
#evidence_df = evidence_df.sample(frac=0.1)

####Filter based on number of observations
badscores = evidence_df[evidence_df['Weighting'] < 1].index
evidence_df.drop(badscores, inplace=True)
evidence_df.index = range(len(evidence_df))


#########################RANDOM FOREST REGRESSOR START##################################


#Make empty dictionaries
RFR_mse_test = dict()
RFR_cv_accuracy = dict()



######ADDING NEW COLUMNS FOR RF REGRESSOR####
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


####SELECTING AND SPLITTING########

feature_subset_RFR = ['Charge','Length','A.Norm', 'R.Norm',
       'N.Norm', 'D.Norm', 'C.Norm', 'Q.Norm', 'E.Norm', 'G.Norm', 'H.Norm',
       'I.Norm', 'L.Norm', 'K.Norm', 'M.Norm', 'F.Norm', 'S.Norm', 'T.Norm',
       'W.Norm', 'Y.Norm', 'V.Norm', 'P.Norm', 'pI', 'Mass', 'pyCharge']
target_value = ['AvgMaxCV']

#Make feature list X that will include the one-hot encoded sequence
X = np.concatenate((evidence_df[feature_subset_RFR], hotarray), axis = 1)
#X = evidence_df[feature_subset1]
y = evidence_df[target_value]

# NOTE: Everything above this could probably be in a separate file

#split the data with 80% in training set
X1, X2, y1, y2 = train_test_split(X, y, random_state=0, train_size=0.80)

y1 = y1.values.ravel() #Flatten the vector so that it is formatted correctly
y2 = y2.values.ravel()

### run single parameter set (using stdin) ###

params_index = int(sys.argv[1])
params_df = pd.read_csv("params.tsv", sep = "\t")
# grab the parameter set based on the index
params = params_df.iloc[params_index]
# initial the dataframe to record model parameters and performance
score_df = pd.DataFrame([params], index=[params_index])
model = RandomForestRegressor(
    max_depth=params['max_depth'],
    min_samples_leaf=params['min_samples_leaf'],
    min_samples_split=params['min_samples_split'],
    n_estimators=params['n_estimators']
)
# get the cross validation scores, keep track of time
start = timeit.default_timer()
scores = cross_val_score(model, X1, y1, cv=5)
stop = timeit.default_timer()
# add the mean
score_df['mean_score'] = round(np.mean(scores),3)
# get the runtime
run_time = stop - start
score_df['run_time'] = round(run_time,3)

# write out the file
outfile_name = "{}.tsv".format(params_index)
score_df.to_csv(outfile_name, index=False, sep = "\t") # index already built into table

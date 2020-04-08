
########Here import all of package pieces I need
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

#split the data with 80% in training set
X1, X2, y1, y2 = train_test_split(X, y, random_state=0, train_size=0.80)

y1 = y1.values.ravel() #Flatten the vector so that it is formatted correctly
y2 = y2.values.ravel()


###WORKING THROUGH PARAMETER GRID

#Had to try to shorten parameters, so these are what I ended up with (6/10/19):
gscRFR = GridSearchCV(
    estimator=RandomForestRegressor(),
    param_grid={'max_depth': (10, 20, 40),
            'n_estimators': (200, 500, 1000),
            'min_samples_leaf': (5, 10, 30),
            'min_samples_split': (5, 50, 200),
            'max_features': (0.25, 0.5, 0.75)},
    cv=5, scoring='neg_mean_squared_error', verbose=3)
#If had unlimited time would run these parameters
        #'max_depth': (None, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150),
        #'n_estimators': (100, 500, 1000, 4000, 8000),
        #'min_samples_leaf': (1, 5, 10, 20, 40, 100, 200, 500, 700, 1000),
        #'min_samples_split': (2, 5, 10, 20, 30, 40, 50, 100, 150, 300, 500, 1000),
        #'max_features': ("auto","log2", "sqrt", 0.01, 0.1, 0.25, 0.5, 0.75, 0.9)},




#Pull out grid search results into a data frame
RFR_gsc_results = gscRFR.fit(X1,y1)
RFR_cvresults_frame = pd.DataFrame(RFR_gsc_results.cv_results_)
#RFR_mse_both_frame = RFR_cvresults_frame[['params', 'mean_train_score', 'mean_train_score']]



#for row in range(0, RFR_cvresults_frame.shape[1]):

    #myforest = RandomForestRegressor(max_depth = RFR_cvresults_frame['param_max_depth'][row], n_estimators = RFR_cvresults_frame['param_n_estimators'][row],
    #                                min_samples_leaf= RFR_cvresults_frame['param_min_samples_leaf'][row], min_samples_split = RFR_cvresults_frame['params_min_samples_split'][row],
    #                                max_features = RFR_cvresults_frame['param_max_features'][row])
    #myforest.fit(X1,y1)
    #y2_rfmodel = myforest.predict(X2)
    #Store metrics in dictionary
    #RFR_mse_test["max_depth: " + str(RFR_cvresults_frame['param_max_depth'][row]) + "; n_estimators: " +  str(RFR_cvresults_frame['param_n_estimators'][row]) +
    #     "; min_samples_leaf: " + str(RFR_cvresults_frame['param_min_samples_leaf'][row]) + "; min_samples_split: " + str(RFR_cvresults_frame['params_min_samples_split'][row]) +
    #     "; max_features: " +  str(RFR_cvresults_frame['param_max_features'][row])] = np.sqrt(mean_squared_error(y2, y2_rfmodel))
    #scores = cross_val_score(myforest, X1, y1, cv=5)    #Does it matter whether this looks at X2 or X1?
    #RFR_cv_accuracy["max_depth: " + str(RFR_cvresults_frame['param_max_depth'][row]) + "; n_estimators: " +  str(RFR_cvresults_frame['param_n_estimators'][row]) +
    #     "; min_samples_leaf: " + str(RFR_cvresults_frame['param_min_samples_leaf'][row]) + "; min_samples_split: " + str(RFR_cvresults_frame['params_min_samples_split'][row]) +
    #     "; max_features: " +  str(RFR_cvresults_frame['param_max_features'][row])] = "Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2)



#Turn dictionaries into dataframes where the keys are the row names
#RFR_mse_test_frame = pd.DataFrame.from_dict(RFR_mse_test, orient = index)
#RFR_cv_accuracy_frame = pd.DataFrame.from_dict(RFR_cv_accuracy, orient = index)

#########################RANDOM FOREST REGRESSOR END##################################

pd.DataFrame(RFR_cvresults_frame).to_csv("RFR_gridserach_results.csv")
#pd.DataFrame(RFR_mse_test_frame).to_csv("RFR_MSE_test.csv")
#pd.DataFrame(RFR_cv_accuracy_frame).to_csv("RFR_CV_Accuracy_training.csv")

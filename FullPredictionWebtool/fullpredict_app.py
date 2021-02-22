## So should be able to just load models and make predictions from sequence

from flask import Flask, render_template, request
import numpy as np
import pandas as pd
import sklearn
import tensorflow.keras.backend as K
from tensorflow import keras
import tensorflow as tf
from numpy import mean
from numpy import std
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
import pickle
from tensorflow.keras.metrics import binary_accuracy
from sklearn.metrics import roc_auc_score, fbeta_score, recall_score, precision_score, accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import roc_auc_score
from fastai.metrics import accuracy_thresh
from fastai.metrics import fbeta
from fastai.torch_core import np2model_tensor
from sklearn.metrics import accuracy_score

import re

#Functions for adding features including pI and one-hot encoding
#BEWARE, VOCAB LENGTH IS FIXED RATHER THAN BASING OFF LONGEST SEQUENCE IN INPUT FOR ONE-HOT
from pyteomics import mass
from pyteomics import parser
from pyteomics import electrochem
import numpy as np
from numpy import array
from numpy import argmax
from keras.utils import to_categorical
from keras_preprocessing.sequence import pad_sequences

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

def addfeatures(featurestable, seqlabel = 'Sequence'):
    Mass = list()
    pI = list()
    Charge = list()

    for i in range(0, featurestable.shape[0]):
        ps = parser.parse(featurestable[seqlabel][i], show_unmodified_termini=True)

        Mass.append(mass.calculate_mass(parsed_sequence=ps))
        Charge.append(electrochem.charge(ps, 2.5))
        pI.append(electrochem.pI(featurestable[seqlabel][i]))


    featurestable['pyMass'] = Mass
    featurestable['pI'] = pI
    featurestable['pyCharge'] = Charge

    return(featurestable)


def simpleOneHot(data_frame, sequenceTag = 'ModSequence', alphabet = 'ACDEFGHIKLMNPQRSTVWY'):
    #Start by finding the max and calculating needed vector length
    VEC_LENGTH = 50 * len(alphabet) #max(data_frame['Length']) * len(alphabet)

    #Define what residues are possible
    AMINO_ACIDS = alphabet

    #TURNING CHARACTERS INTO INTEGERS
    # Map character keys to integer values in a dictionary, then map integer keys to character values to revers transform
    char_to_int = dict((c, i) for i, c in enumerate(AMINO_ACIDS))   #character keys to integer values
    int_to_char = dict((i, c) for i, c in enumerate(AMINO_ACIDS))   #integer keys to character values


    hotlist = list()
    #Build out the rest of the sequences' one-hot arrays

    for i in range(0, data_frame.shape[0]):

        pep = data_frame[sequenceTag][i]
        #print(pep)
        integer_encode = [char_to_int[char] for char in pep]
        encoded = to_categorical(integer_encode, num_classes=22)
        flatencode = encoded.flatten()

        #numzeros = VEC_LENGTH - len(flatencode)
        #flatencode = np.append(flatencode, [[0] * numzeros])

        hotlist.append(flatencode)

    padded = pad_sequences(hotlist, padding= 'post', maxlen=VEC_LENGTH)

    hotarray = np.array(padded)

    hotarray.shape
    return(hotarray)

def fbeta2(y_true, y_pred, threshold_shift=0):
    beta = 2

    # just in case of hipster activation at the final layer
    y_pred = K.clip(y_pred, 0, 1)

    # shifting the prediction threshold from .5 if needed
    y_pred_bin = K.round(y_pred + threshold_shift)

    tp = K.sum(K.round(y_true * y_pred_bin)) + K.epsilon()
    fp = K.sum(K.round(K.clip(y_pred_bin - y_true, 0, 1)))
    fn = K.sum(K.round(K.clip(y_true - y_pred, 0, 1)))

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    beta_squared = beta ** 2
    return (beta_squared + 1) * (precision * recall) / (beta_squared * precision + recall + K.epsilon())

def f2(y_true, y_pred):
    y_pred = K.round(y_pred)
    tp = K.sum(K.cast(y_true*y_pred, 'float'), axis=0)
    # tn = K.sum(K.cast((1-y_true)*(1-y_pred), 'float'), axis=0)
    fp = K.sum(K.cast((1-y_true)*y_pred, 'float'), axis=0)
    fn = K.sum(K.cast(y_true*(1-y_pred), 'float'), axis=0)

    p = tp / (tp + fp + K.epsilon())
    r = tp / (tp + fn + K.epsilon())

    f1 = 5*p*r / (4*p+r+K.epsilon())
    f1 = tf.where(tf.math.is_nan(f1), tf.zeros_like(f1), f1)
    return K.mean(f1)



#app = Flask('Prediction_Webtool')

#@app.route('/')
#def show_predict_stock_form():
#    return render_template('predictorform.html')

#@app.route('/result', methods = ['POST'])


#write your function that loads the model
#Load RandomForestClassifier
loaded_rfc = pickle.load(open('/Users/justinmcketney/Desktop/FAIMS_MachineLearning/standard_models/RFC_TrainedOnFullHuman.sav', 'rb'))
#Load Neuaral Network
mload = keras.models.load_model('/Users/justinmcketney/Desktop/FAIMS_MachineLearning/neural_net/final_all_human_model.h5', custom_objects={'f2':f2, 'fbeta2':fbeta2})


def generate_predictions(input_seq):
    #GENERATING FULL FEATURES FOR RANDOM FOREST
    ### Generate charge states, just do all for right now, we'll filter later
    input_seq_df = pd.DataFrame({'Charge' : [2,3,4,5], 'Sequence' : [input_seq] * 4})
    input_seq_df['SeqCharge'] =  input_seq_df['Charge'].astype(str) + input_seq_df['Sequence']
    input_seq_df['Length'] = len(input_seq_df['Sequence'])
    input_seq_df['ModSequence'] = input_seq_df['Sequence']
    #### Adding other features from pyteomics
    input_seq_df = addfeatures(input_seq_df)
    #### Creating one-hot array for our peptides (this might be problem because of dimensions)
    input_seq_hotarray = simpleOneHot(data_frame=input_seq_df, alphabet='ACDEFGHIKLMNPQRSTVWYam')
    #### Creating final feature array for RFC
    feature_subset = ['Charge', 'Length', 'pyMass', 'pI']
    RFC_x = np.concatenate((input_seq_df[feature_subset], input_seq_hotarray), axis = 1)

    #GENERATING FULL FEATURES FOR NEURAL NET
    ### combine all the letters into a long string, take the set to find the unique values, add 'END' (for use with one-hot), then get length
    seq = input_seq_df['SeqCharge']
    vocab = set(''.join([str(i) for i in seq]))
    vocab.add('END')
    len_vocab = len(vocab)
    ### char index is static
    char_index = {'2': 0, '3': 1, 'F': 2, 'a': 3, 'E': 4, 'T': 5, 'M': 6, '5': 7, 'm': 8, 'R': 9, 'END': 10, 'V': 11, 'A': 12, 'K': 13, 'I': 14, 'G': 15, 'W': 16, 'P': 17, 'Q': 18, 'D': 19, '4': 20, 'C': 21, 'N': 22, 'L': 23, 'S': 24, 'Y': 25, 'H': 26}
    ### BEWARE fixing max length at 51
    ### maxlen = max([len(x) for x in df.SeqCharge])
    maxlen = 51
    NN_x = []
    x_name = [str(i)[0:maxlen] for i in seq]
    for i in x_name:
        tmp = [char_index[j] for j in str(i)]
        for k in range(0,maxlen - len(str(i))):
            tmp.append(char_index["END"])
        NN_x.append(tmp)
    NN_x = np.asarray(NN_x)

    #GENERATING PREDICTIONS AND COMBINING INTO DATAFRAME
    CV_settings = ['CV20', 'CV25', 'CV30', 'CV35','CV40', 'CV45', 'CV50', 'CV55','CV60', 'CV65', 'CV70', 'CV75','CV80', 'CV85', 'CV90', 'CV95']
    nn_predictions_output = mload.predict(NN_x)
    rfc_predictions_output = loaded_rfc.predict_proba(RFC_x)
    ensemble_predictions = (rfc_predictions_output + nn_predictions_output)/ 2.0

    ## Make three data frames and then concatenate
    ### NN dataframe
    nn_predict_df = pd.DataFrame(nn_predictions_output, columns= CV_settings)
    nn_predict_df['Peptide Ion'] = input_seq_df['SeqCharge']
    nn_predict_df['Model'] = 'Neural Net'
    ### RFC dataframe
    rfc_predict_df = pd.DataFrame(rfc_predictions_output, columns= CV_settings)
    rfc_predict_df['Peptide Ion'] = input_seq_df['SeqCharge']
    rfc_predict_df['Model'] = 'Random Forest'
    ### Ensemble dataframe
    ens_predict_df = pd.DataFrame(ensemble_predictions, columns= CV_settings)
    ens_predict_df['Peptide Ion'] = input_seq_df['SeqCharge']
    ens_predict_df['Model'] = 'Ensemble'

    concatenated_predictions = pd.concat([rfc_predict_df, nn_predict_df, ens_predict_df], axis=0)

    ## Sort the columns and rows
    cols = concatenated_predictions.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    concatenated_predictions = concatenated_predictions[cols]
    prediction_clf = concatenated_predictions.sort_values(by = ['Peptide Ion', 'Model'])

    return prediction_clf

def generate_table(dataframe, max_rows=10):
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in dataframe.columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
            ]) for i in range(min(len(dataframe), max_rows))
        ])
    ])

#user_seq = "PEPTIDE"
df = generate_predictions("PEPTIDE")

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1(children = 'Peptide Submission'),
    dcc.Input(id='user_seq', type='text', value='PEPTIDE'),
    html.Button(id='submit-button-state', n_clicks=0, children='Submit'),
    html.H1(children = 'Prediction Probabilities'),
    html.Div(id ='output-state')
    ])

@app.callback(Output('output-state', 'children'),
              Input('submit-button-state', 'n_clicks'),
              State('user_seq', 'value'))

def update_output(n_clicks, user_seq):

    #Convert the user_seq to upper case
    user_seq_upper = user_seq.upper()
    #Make list of unacceptable characters
    regex_reject = re.compile('[@_!#$%^&*()<>?/\|}{~:]|[0-9]|[BJOUXZ]')
    if(regex_reject.search(user_seq_upper) != None):
        new_fig = 'There is an non amino acid character in your sequence'
    #Check if sequence is too short
    if(len(user_seq_upper) < 7):
        new_fig = 'The sequence supplied is too short for reliable prediction (< 7 residues)'

    #Check if sequence is too long
    if(len(user_seq_upper) > 48):
        new_fig = 'The sequence supplied is too long for reliable prediction (max 48 residues)'

    #If everything is just right (Goldilox style) generate the predicitons
    # and round to two decimal places
    if(regex_reject.search(user_seq_upper) == None and len(user_seq_upper) > 6 and len(user_seq_upper) < 49):
        new_df = generate_predictions(user_seq_upper)
        round_df = new_df.round(3)
        new_fig = generate_table(round_df)

    return new_fig

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')

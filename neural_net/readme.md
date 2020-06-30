# minimum working example of neural network training

Currently just uses embedding layer and dual LSTM - binaryAccuracy of 0.9 ish on the validation set

##### BEST HYPERPARAMETERS:

embedding neurons = 96
LSTM neurons = 96
Dropout1 = 0.32255
Dense layer neurons = 64
dropout2 = 0.02152
learning rate = 0.001758
best # training epochs = 52


##### 30% human test set performance:

{'loss': 0.19521527745990358,
 'acc': 0.9072674,
 'prec': 0.66095227,
 'recall': 0.5580575,
 'roc_auc': 0.9408676,
 'f2': 0.57525325}

make sure to use custom_objects= {'f2': f2, 'fbeta2':fbeta2} when importing the model or it will throw an error!

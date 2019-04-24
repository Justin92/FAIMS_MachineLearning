
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# Load data
norm_data_path = "/Volumes/GENERAL/File_exchange/Ian/From_Justin/NormCountsFeatures_MaxCV.csv"
evidence_df = pd.read_csv(norm_data_path, low_memory=False)

# select and decompose normalized count data
normalized_aa_counts = ['A.Norm', 'R.Norm',
       'N.Norm', 'D.Norm', 'C.Norm', 'Q.Norm', 'E.Norm', 'G.Norm', 'H.Norm',
       'I.Norm', 'L.Norm', 'K.Norm', 'M.Norm', 'F.Norm', 'S.Norm', 'T.Norm',
       'W.Norm', 'Y.Norm', 'V.Norm', 'U.Norm', 'O.Norm', 'P.Norm']

normalized_count_df = evidence_df[normalized_aa_counts]

# reduce dimensionality or normalized AA counts with PCA
pca = PCA(n_components=10)
decomposed_df = pca.fit_transform(normalized_count_df)

# concat decomposed values with length, charge
recombined_df = pd.concat([evidence_df[['Max Intensity CV', 'Length.x', 'Charge']]
    , pd.DataFrame(decomposed_df]), axis=1)

# set target value and create training data
target_value = ['Max Intensity CV']
# exclude target value (first column) from prepared feature matrix
X = recombined_df[recombined_df.columns[1:]]
y = recombined_df[target_value]

#split the data with 50% in each set
X1, X2, y1, y2 = train_test_split(X, y, random_state=0, train_size=0.8)

#Flatten the vector so that it is formatted correctly
y1 = y1.values.ravel()
y2 = y2.values.ravel()

# set hyperparameters and fit model
forest = RandomForestRegressor(200)
forest.fit(X1, y1)

# make predictions based on fitted model
y2_model = forest.predict(X2)
predicted = cross_val_predict(forest1, X1, y1, cv=5, n_jobs=-1)

# evaluate performance with RMSE
np.sqrt(mean_squared_error(y2, y2_model))

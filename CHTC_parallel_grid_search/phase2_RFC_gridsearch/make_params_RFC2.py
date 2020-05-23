from sklearn.model_selection import cross_val_predict, GridSearchCV
from sklearn.model_selection import ParameterGrid
from sklearn.ensemble import RandomForestClassifier
import pandas as pd

#Code here defines series of paramter dictionaries that in the GridSearchCV
#Grid search setting essentially just dictionary of dictionaries
gscRFC = GridSearchCV(
    estimator=RandomForestClassifier(),
    param_grid={
        'max_depth': (30, 60, 90, 120), #15, 20, 25, None), #, 30, 40, 50, 60, 80, 100, 150),
        'n_estimators': (50, 100, 200, 300), #500, 1000), #, 4000, 8000),
        'min_samples_split': (5, 10, 20, 40, 60), #10, 20, 30, 40),#, 50, 100, 150, 300, 500, 1000),
        'max_features': ( 0.2, 0.4, 0.5, 0.75, 0.99)}, #0.25, 0.5, 0.75, 0.9)},

    cv=5, scoring='neg_mean_squared_error', verbose=3)

# produce a list of param dictionaries because gscRFR.param_grid is a dictionary
#ParameterGrid makes it into a ParameterGrid type object
param_list = list(ParameterGrid(gscRFC.param_grid))

# save params as tsv
#Do this by making it first a ParameterGrid, then list then dataframe
params_df = pd.DataFrame(list(ParameterGrid(gscRFC.param_grid)))
params_df.index.name = "index"
params_df.to_csv("RFC_params2.tsv", sep = "\t")
params_indices = pd.DataFrame(params_df.index)
params_indices.to_csv("RFC_params_indices2", sep = "\t", index = False)

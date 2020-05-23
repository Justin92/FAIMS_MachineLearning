from sklearn.model_selection import cross_val_predict, GridSearchCV
from sklearn.model_selection import ParameterGrid
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
import pandas as pd



#Code here defines series of paramter dictionaries that in the GridSearchCV
#Grid search setting essentially just dictionary of dictionaries
gscRFC = GridSearchCV(
    estimator=XGBClassifier(),
    param_grid={
        'max_depth': (2, 3),
        #'booster': ('dart', 'gbtree'),
        #'learning_rate': (0.1, 0.3, 0.5),
        'max_delta_step': (0, 1),
        'scale_pos_weight': (0.01, 0.05),
        #'subsample': (0.75), #, 1.0),
        #'tree_method': ("auto", "hist")
        },

    cv=5, scoring='neg_mean_squared_error', verbose=3)

# produce a list of param dictionaries because gscRFR.param_grid is a dictionary
#ParameterGrid makes it into a ParameterGrid type object
param_list = list(ParameterGrid(gscRFC.param_grid))

# save params as tsv
#Do this by making it first a ParameterGrid, then list then dataframe
params_df = pd.DataFrame(list(ParameterGrid(gscRFC.param_grid)))
params_df.index.name = "index"
params_df.to_csv("XGBC_params.tsv", sep = "\t")
params_indices = pd.DataFrame(params_df.index)
params_indices.to_csv("XGBC_params_indices", sep = "\t", index = False)

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import sys, os

import numpy as np
import pandas as pd
from scipy.special import comb
from scipy import stats
import scipy.cluster.hierarchy as hac
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import metrics
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import StackingClassifier

from sklearn.metrics import balanced_accuracy_score
from joblib import dump, load
from sklearn.preprocessing import LabelEncoder


# https://www.kaggle.com/code/rafjaa/dealing-with-very-small-datasets#t1
# 
# https://medium.com/rants-on-machine-learning/what-to-do-with-small-data-d253254d1a89




# only use the first 90 days because that's what we have from Yara
#timepoints_10cc_transfer0 = refdf_10cc.loc[
#    refdf_10cc.transfer.isin([0]), # & refdf_10cc.day.le(91), 
#    'day'].unique()
#timepoints_10cc_transfer0


# # interpolation to make all groups have the same timepoints

def interpolate_series(x, interpolate_cols, index_col, timepoints,):
    x1 = x[[index_col] + interpolate_cols].copy()
    x1.set_index(index_col, inplace=True)
    x1 = x1.reindex(x1.index.union(timepoints))
    x1 = x1.interpolate(method='index',limit_direction='both')
    x1 = x1.reindex(timepoints)    
    return x1



def interpolate_dataframe(df, groupby_cols, interpolate_cols, index_col, timepoints,):
    func = lambda x: interpolate_series(x, interpolate_cols, index_col, timepoints,)

    refdf_interpolate = df.groupby(groupby_cols).apply(func)
    refdf_interpolate.reset_index(inplace=True)
    return refdf_interpolate

def df2clipped(df, X_cols=['ref_Bp[N]', 'ref_Bp[C]'], max_X_day=91, clip_value=1):
    # only use the first 90 days because that's what we have from Yara
    df1 = df.loc[df.day.le(max_X_day)].copy()
    for i in X_cols:
        df1[i] = df1[i].clip(lower=clip_value)
    return df1

def df2X(df, groupby_cols, X_cols=['ref_Bp[N]', 'ref_Bp[C]'], index_col='day', max_X_day=91, clip_value=1):
    # only use the first 90 days because that's what we have from Yara
    X = df.pivot_table(index=groupby_cols, values=X_cols,columns=index_col)
    return X

def X2logged(X):
    X_logged = np.log(X)
    return X_logged

#    y = X.index.get_level_values(0)



# # Compute features



def _compute_max_features(df, groupby_cols, nutrient, biomass_prefix, index_col='day'):
    bcol=f'{biomass_prefix}[{nutrient}]'
    refmaxbp_df = df.loc[df.groupby(groupby_cols)[bcol].idxmax()].reset_index(drop=True)
    refmaxbp_df = refmaxbp_df.set_index(groupby_cols)
    refmaxbp_df = refmaxbp_df[[index_col, bcol]]
    refmaxbp_df.rename(columns={
        index_col : f'max_{nutrient}day', 
        bcol  : f'max_{nutrient}biomass',
    }, inplace=True)
    return refmaxbp_df

def _compute_mean_features(df, groupby_cols, nutrient, biomass_prefix):
    bcol=f'{biomass_prefix}[{nutrient}]'
    min_day = 30
    max_day = 60
    lterm_df = df.loc[df.day.ge(min_day) & df.day.le(max_day)]
    
    reflterm_df = lterm_df.groupby(groupby_cols)[bcol].agg(['mean', 'std', 'median'])
    reflterm_df.rename(columns={
        'mean' : f'mean_{nutrient}biomass',
        'median' : f'median_{nutrient}biomass',
        'std' : f'std_{nutrient}biomass',
    }, inplace=True)
    return reflterm_df

def _compute_lastday_features(df, groupby_cols, nutrient, biomass_prefix):
    bcol=f'{biomass_prefix}[{nutrient}]'
    lod_lastday_threshold = 2
    refmaxday_df = df.loc[df[bcol].ge(lod_lastday_threshold)].groupby(groupby_cols).day.max()
    #refmaxday_df.rename(columns=dict(day=f'last_day{nutrient}'), inplace=True)
    refmaxday_df.name = f'last_day{nutrient}'
    return refmaxday_df


def compute_features(df, groupby_cols, biomass_prefix='ref_Bp'):
    ''' compute features for random forest
    '''
    # need to do that so idxmax will give the correct results
    df1 = df.reset_index(drop=True)

    df_list = (
        [_compute_max_features(df1, groupby_cols, nutrient, biomass_prefix) for nutrient in 'NC'] +
        [_compute_mean_features(df1, groupby_cols, nutrient, biomass_prefix) for nutrient in 'NC'] +
        [_compute_lastday_features(df1, groupby_cols, nutrient, biomass_prefix) for nutrient in 'NC']
    )
    df_merge = df_list[0].join(df_list[1:])
        
    df_merge.fillna(0, inplace=True) # for last day
    return df_merge
    


# log the biomass
def add_log_cols(df, biomass_prefix='ref_Bp'):
    lod_threshold = 1
    for nutrient in 'NC':
        df[f'log_{nutrient}biomass'] = np.log(df[f'{biomass_prefix}[{nutrient}]'].clip(lower=lod_threshold))
    return df



def _X_smt_to_features(X_smt):
    x = X_smt.T.melt(ignore_index=False, var_name='smt_id',).reset_index()
    x = x.pivot(index=['smt_id', 'day'], values='value', columns='level_0').reset_index()
    return _compute_features(x, groupby_col='smt_id', biomass_prefix='ref_Bp')




def df2finalX(df, groupby_cols, log_mode=True):
    dfclipped =  df2clipped(df)
    Xclipped = df2X(dfclipped, groupby_cols)
    if log_mode:
        X_logged = X2logged(Xclipped)
    else:
        X_logged = Xclipped

    feat_df = compute_features(dfclipped, groupby_cols)

    X_logged.columns = [f'{col}_{day:2.1f}' for col,day in X_logged.columns.values]
    X_final = feat_df.join(X_logged)

    logistic_Nfeatures = [c for c in X_final.columns if c.startswith('ref_Bp[N]')]
    logistic_Cfeatures = [c for c in X_final.columns if c.startswith('ref_Bp[C]')]
    forest_features = feat_df.columns
    return (X_final, forest_features, logistic_Nfeatures, logistic_Cfeatures)
    

#def X_logged_to_final_X(X_logged):
#def df2finalX(df, groupbycols):
# sim_df = sim_df[['run_id', 'day', 'Bptotal[N]','Bptotal[C]']]
# sim_df.rename(columns={'Bptotal[N]':'ref_Bp[N]','Bptotal[C]':'ref_Bp[C]'}, inplace=True)
# for c in ['ref_Bp[N]', 'ref_Bp[C]']:
    # sim_df[c] = np.log10(sim_df[c].clip(lower=1))
# sim_df_filtered = sim_df.loc[sim_df.day.le(91)]
# X_sim_logistic = sim_df_filtered.pivot_table(index='run_id', values=['ref_Bp[N]', 'ref_Bp[C]'],columns='day')
# X_sim_logistic1.loc[X_sim_logistic1.isna().sum(axis=1).ge(1)]
# # meed to do that so idxmax will give the correct results
# sim_df_filtered = sim_df_filtered.reset_index(drop=True)
# sim_groupby_col = ['run_id']
# psim_feature_df = _compute_features(sim_df_filtered, sim_groupby_col, 'ref_Bp')
# X_sim_logistic1.columns = [f'{col}_{day:2.1f}' for col,day in X_sim_logistic1.columns.values]
# X_sim = psim_feature_df.join(X_sim_logistic1)



  
# # Stacking classifier
# Classifiers
# The goal is to create multiple classifiers:
# 1. randomforest based on generated features: last day, max day, max value, mean late biomass (N & C), std late biomass (N & C)
# 2. logical regression based on N biomass
# 3. logical regression based on C biomass
# 4. tslearn based on N biomass 
# 5. tslearn based on C biomass 
# 
# and than combine the results using a stacking classifier

    


# Classifiers
def build_classifier(forest_features, logistic_Nfeatures, logistic_Cfeatures, y_train):
    le = LabelEncoder()
    le.fit_transform(y_train), le.classes_

    class_weight_map = {'Axenic' : 100, 'Inhibited' : 10000, 'Other': 1, 'Strong' : 10000, 'Sustained' : 10000, 'Weak': 100}
    # TODO le.classes_
    class_weight = {i : class_weight_map[c] for i,c in enumerate(le.classes_)}

    clf_features = RandomForestClassifier(ccp_alpha=0.05, class_weight=class_weight)
    clf_logisticN = LogisticRegression(penalty='l2', C=0.1, max_iter=10000, class_weight=class_weight,)
    clf_logisticC = LogisticRegression(penalty='l2', C=0.1, max_iter=10000, class_weight=class_weight,)
    clf_meta = LogisticRegression()

    pipe_logisticN = Pipeline([
        ('select', ColumnTransformer([('sel', 'passthrough', logistic_Nfeatures)], remainder='drop')),  
        ('scale', StandardScaler()),    
        ('clf', clf_logisticN)
    ],
        memory='/tmp/Osnat/sklearn_cache',
    )
    pipe_logisticC = Pipeline([
        ('select', ColumnTransformer([('sel', 'passthrough', logistic_Cfeatures)], remainder='drop')),
        ('scale', StandardScaler()),    
        ('clf', clf_logisticC)
    ],
        memory='/tmp/Osnat/sklearn_cache',
    )
    pipe_forest = Pipeline([
        ('select', ColumnTransformer([('sel', 'passthrough', forest_features)], remainder='drop')),  
        ('scale', StandardScaler()),    
        ('clf', clf_features)
    ],
        memory='/tmp/Osnat/sklearn_cache',
    )



    stack = StackingClassifier(
        estimators=[
            ('logisticN', pipe_logisticN),
            ('logisticC', pipe_logisticC),
            ('ForestFeatures', pipe_forest),
        ],
        final_estimator=clf_meta,
        n_jobs = -1, verbose=True
    )
    return stack






# # classify simulations



# sim_df = sim_df[['run_id', 'day', 'Bptotal[N]','Bptotal[C]']]
# sim_df.rename(columns={'Bptotal[N]':'ref_Bp[N]','Bptotal[C]':'ref_Bp[C]'}, inplace=True)




# for c in ['ref_Bp[N]', 'ref_Bp[C]']:
    # sim_df[c] = np.log10(sim_df[c].clip(lower=1))


# sim_df_filtered = sim_df.loc[sim_df.day.le(91)]

# X_sim_logistic = sim_df_filtered.pivot_table(index='run_id', values=['ref_Bp[N]', 'ref_Bp[C]'],columns='day')


# X_sim_logistic1.loc[X_sim_logistic1.isna().sum(axis=1).ge(1)]

# # meed to do that so idxmax will give the correct results
# sim_df_filtered = sim_df_filtered.reset_index(drop=True)

# sim_groupby_col = ['run_id']

# psim_feature_df = _compute_features(sim_df_filtered, sim_groupby_col, 'ref_Bp')

# X_sim_logistic1.columns = [f'{col}_{day:2.1f}' for col,day in X_sim_logistic1.columns.values]

# X_sim = psim_feature_df.join(X_sim_logistic1)

def save_ML_model(stack, fpath = '10CC_ML_classifier.joblib'):
    dump(stack, fpath)


def load_ML_model(fpath = '10CC_ML_classifier.joblib'):
    stack = load(fpath)
    return stack
    
def ML_model_predict_proba(stack, X_sim):
    y_sim_pred_prop = stack.predict_proba(X_sim)
    y_sim_pred = stack.predict(X_sim)
    max_sim_prob = np.amax(y_sim_pred_prop, axis=1)
    # df_predicted_prob = pd.DataFrame(y_sim_pred_prop, columns=stack.classes_)
    # df_predicted_prob['run_id'] =  X_sim.index
    # df_predicted_prob['y_pred'] =  y_sim_pred
    # df_predicted_prob['max_prob'] =  max_sim_prob
    # mdf_predicted_prob = df_predicted_prob.melt(
        # id_vars=['run_id','y_pred', 'max_prob'], 
        # value_vars=stack.classes_, 
        # value_name='prob', var_name='Group')
    df_predicted_classes = pd.DataFrame({
        # 'run_id' : X_sim.index,
        'y_pred' : y_sim_pred,
        'max_prob' : max_sim_prob,
    }, index=X_sim.index,)
    #df_predicted_classes[['idx', 'media', 'which', 'model', 'hash']] = df_predicted_classes.run_id.str.rsplit('_', n=4, expand=True)
    return df_predicted_classes



# TODO
#df_predicted_classes.to_csv('monte_predicted_classes.csv.gz', index=False)


def get_features_dataframe(X_sim, forest_features):
    df_sim_maxday  = X_sim[forest_features].reset_index()
    # df_sim_maxday = df_sim_maxday.melt(id_vars='run_id')
    # df_sim_maxday = pd.merge(df_sim_maxday, df_predicted_classes, on='run_id')
    return df_sim_maxday



#ref_groupby_col = ['train_Group', 'Sample', 'id', 'full name', 'Group', 'Experiment',]



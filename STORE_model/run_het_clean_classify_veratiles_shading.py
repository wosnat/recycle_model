#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import re
import os
import ML_model.run_ML_classification as ml

from scipy.special import comb
from scipy import stats
import scipy.cluster.hierarchy as hac
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import statsmodels.api as sm

from joblib import dump, load



#from model_equations_separate_NC import *
from model_equations_separate_NC_store_numba import *

# represent FL 0.1
NBIOMASS_LOD_NOT_GROWING = 2.34502821

# pro99_mode = False 
# which_organism = 'ponly'
# organism_to_tune = 'PRO'
# (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
# ) = get_constants_per_organism(pro99_mode, which_organism)
# for model in ['OVERFLOW', 'MIXOTROPH', 'ROS', 'MIN', 'EXOENZYME']:
    # params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)
    


# cleanup
def cleanup_session(session_id, sim_df, mse_df, sum_df, out_dpath, disable_C2N):
    pro99_mode = False 
    which_organism = 'all'
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
    ) = get_constants_per_organism(pro99_mode, which_organism)

    dpath = 'results'


    min_values = sim_df[var_names+['run_id']].groupby('run_id').min().min(axis=1)
    bad_ids_negative_values = set(min_values[min_values< -1e-9].index)
    print('bad_ids_negative_values: ',len(bad_ids_negative_values))

    max_values = sim_df[var_names+['run_id']].groupby('run_id').max().max(axis=1)
    bad_ids_toolarge_values = set(max_values[max_values > 1e9].index)
    print('bad_ids_toolarge_values: ',len(bad_ids_toolarge_values))

    bad_ids_missing_compare_points =  mse_df.loc[mse_df.compare_points != 74, 'run_id']
    print('bad_ids_missing_compare_points: ',len(bad_ids_missing_compare_points))

    if disable_C2N:
        problematic_C2N_runids = set()
    else:
        problematic_C2N_runids = sim_df.loc[
            sim_df['Bhtotal[C]'].ge(5) & 
            sim_df['QCh'].le(3.5)& 
            sim_df['day'].ge(10)
        ].run_id.unique()
    print('problematic_C2N_runids: ',len(problematic_C2N_runids))
    
    

    bad_ids = set(bad_ids_missing_compare_points) | set(bad_ids_negative_values) | set(bad_ids_toolarge_values) | set(problematic_C2N_runids)
    print('bad_ids: ',len(bad_ids))

        
    sim_df1 = sim_df.loc[~sim_df['run_id'].isin(bad_ids)].copy()
    mse_df1 = mse_df.loc[~mse_df['run_id'].isin(bad_ids)].copy()
    sum_df1 = sum_df.loc[~sum_df['run_id'].isin(bad_ids)].copy()

    stats_df = pd.DataFrame([{
        'session_id': session_id,
        'Montecarlo samples':mse_df.run_id.nunique(), 
        'problematic_C2N_runids' : len(problematic_C2N_runids),
        'Bad Samples':len(bad_ids), 
        'Montecarlo cleaned':mse_df1.run_id.nunique(), 
    }]) #.to_csv(os.path.join(out_dpath, f'cleanup_stats_{session_id}.csv')
    

    sim_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_df.csv.gz',), index=False)
    mse_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_mse.csv.gz',), index=False)
    sum_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_sum.csv.gz',), index=False)

    return(sim_df1, mse_df1, sum_df1, stats_df)
    

def classify_samples(sim_df1, mse_df):
    model_fpath = os.path.join(os.path.dirname(__file__),'ML_model','shading_10CC_ML_classifier.joblib')
    stack = load(model_fpath) 
    sim_df = sim_df1[['run_id', 'day', 'Bptotal[N]','Bptotal[C]']].copy()
    sim_df.rename(columns={
            'Bptotal[N]':'ref_Bp[N]',
            'Bptotal[C]':'ref_Bp[C]'
        }, inplace=True)


    X_sim, forest_features, logistic_Nfeatures, logistic_Cfeatures = ml.df2finalX(
        sim_df, groupby_cols=['run_id'])

    df_predicted_classes = ml.ML_model_predict_proba(stack, X_sim)
    df_predicted_classes[['idx', 'media', 'which', 'model', 'hash']] = df_predicted_classes.run_id.str.rsplit('_', n=4, expand=True)

    # TODO
    df_predicted_classes['VPRO'] = df_predicted_classes.idx.str.replace(r'.*(vpro.*)_\d+',r'\1', regex=True)

    df_predicted_classes_merged = pd.merge(df_predicted_classes, mse_df, left_on=['run_id', 'y_pred'], right_on=['run_id', 'Group'], how='left')
    df_predicted_classes_merged['RMSE_filled'] = df_predicted_classes_merged['RMSE'].fillna(0)
    df_predicted_classes_merged =df_predicted_classes_merged.reset_index(drop=True)
    df_predicted_classes_merged_min = df_predicted_classes_merged.loc[
        df_predicted_classes_merged.groupby('run_id')['RMSE_filled'].idxmin()]


    # rename Axenic to neutral
    df_predicted_classes_merged_min.loc[df_predicted_classes_merged_min.y_pred.isin(['Axenic']), 'y_pred'] = 'Neutral'
    
    return df_predicted_classes_merged_min



def find_versatile_vpros(df_predicted_classes):
    vpro_df = df_predicted_classes.pivot_table(
        columns = 'y_pred',
        index=['model', 'VPRO'],
        values='run_id', 
        aggfunc='count', fill_value=0,
    )

    poscolumns = [c for c in ['Strong','Sustained'] if c in vpro_df.columns]
    negcolumns = [c for c in ['Inhibited','Weak'] if c in vpro_df.columns]
    if poscolumns:
        vpro_df['pos_interaction'] = vpro_df[poscolumns].sum(axis=1) 
    else:
        vpro_df['pos_interaction'] = 0

    if negcolumns:
        vpro_df['neg_interaction'] = vpro_df[negcolumns].sum(axis=1)
    else:
        vpro_df['neg_interaction'] = 0
                 
    vpro_df['Versatile'] = vpro_df['pos_interaction'].ge(1) & vpro_df['neg_interaction'].ge(1) 
    vpro_df = vpro_df.reset_index()

    return vpro_df




        

if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='cleanup PONLY runs and produce VPROS.')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--indpath", help="input dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--disable_C2N", help="disable the cleanup C/N filtering",  action="store_true")
    
    
    args = parser.parse_args()
    out_dpath = args.outdpath
    if out_dpath != '':
        os.makedirs(out_dpath, exist_ok=True)

    dpath = args.indpath
    
    session_id = args.run_id
    sim_df = pd.read_csv(os.path.join(dpath,f'{session_id}_df.csv.gz',))
    mse_df = pd.read_csv(os.path.join(dpath,f'{session_id}_mse.csv.gz',))
    sum_df = pd.read_csv(os.path.join(dpath,f'{session_id}_sum.csv.gz',))

    (sim_df, mse_df, sum_df, stats_df) = cleanup_session(session_id, sim_df, mse_df, sum_df, out_dpath, args.disable_C2N)
    
    df_predicted_classes = classify_samples(sim_df, mse_df)
    vpro_df = find_versatile_vpros(df_predicted_classes)
    
    df_predicted_classes.to_csv(os.path.join(out_dpath, f'predicted_classes_{session_id}.csv.gz'), index=False)
    vpro_df.to_csv(os.path.join(out_dpath, f'versatile_vpros_{session_id}.csv'), index=False)

    stats_df.to_csv(os.path.join(out_dpath, f'cleanup_stats_{session_id}.csv'))


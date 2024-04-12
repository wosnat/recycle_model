#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import re
import os

#from model_equations_separate_NC import *
from model_equations_separate_NC_store_numba import *


morder = ['MIN', 'OVERFLOW', 'MIXOTROPH', 'EXOENZYME', 'ROS',]


refdf = pd.read_excel('reference_10cc_axenic.xlsx')
refp99df = pd.read_excel('reference_pro99_axenic.xlsx')



# pro99_mode = False 
# which_organism = 'ponly'
# organism_to_tune = 'PRO'
# (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
# ) = get_constants_per_organism(pro99_mode, which_organism)
# for model in ['OVERFLOW', 'MIXOTROPH', 'ROS', 'MIN', 'EXOENZYME']:
    # params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)
    
session_id = 'ponly_monte_add_ROS_round2'


# cleanup
def cleanup_session(session_id, sim_df, mse_df, sum_df, out_dpath):
    pro99_mode = False 
    which_organism = 'ponly'
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
    ) = get_constants_per_organism(pro99_mode, which_organism)

    dpath = 'results'

    sum_df[['id', 'media', 'which', 'model', 'hash']] = sum_df.run_id.str.rsplit('_', n=4, expand=True)
    mse_df[['id', 'media', 'which', 'model', 'hash']] = mse_df.run_id.str.rsplit('_', n=4, expand=True)
    sim_df[['id', 'media', 'which', 'model', 'hash']] = sim_df.run_id.str.rsplit('_', n=4, expand=True)

    sum_df['idx'] = sum_df['id'] + '_' + sum_df['model']
    mse_df['idx'] = mse_df['id'] + '_' + mse_df['model']
    sim_df['idx'] = sim_df['id'] + '_' + sim_df['model']

    sim_df.drop(columns=['id', 'which', 'hash'], inplace=True)
    mse_df.drop(columns=['id', 'which', 'hash'], inplace=True)
    sum_df.drop(columns=['id', 'which', 'hash'], inplace=True)

    min_values = sim_df[var_names+['idx']].groupby('idx').min().min(axis=1)
    bad_ids_negative_values = set(min_values[min_values< -1e-9].index)
    print('bad_ids_negative_values: ',len(bad_ids_negative_values)

    max_values = sim_df[var_names+['idx']].groupby('idx').max().max(axis=1)
    bad_ids_toolarge_values = set(max_values[max_values > 1e9].index)
    print('bad_ids_toolarge_values: ',len(bad_ids_toolarge_values)


    mse_df['ref_compare_points'] = mse_df.media.map({'lowN' : 74, 'pro99': 58})
    bad_ids_missing_points = set(mse_df.loc[mse_df.ref_compare_points != mse_df.compare_points, 'idx' ])
    print('bad_ids_missing_points: ',len(bad_ids_missing_points)

    bad_ids = set(bad_ids_missing_points) | set(bad_ids_negative_values) | set(bad_ids_toolarge_values)
    print('bad_ids: ',len(bad_ids)

    pd.Dataframe({
        'Montecarlo samples'=mse_df.idx.nunique(), 
        'Bad Samples'=len(bad_ids), 
        
    sim_df1 = sim_df.loc[~sim_df['idx'].isin(bad_ids)].copy()
    mse_df1 = mse_df.loc[~mse_df['idx'].isin(bad_ids)].copy()
    sum_df1 = sum_df.loc[~sum_df['idx'].isin(bad_ids)].copy()

    stats_df = pd.Dataframe({
        'session_id': session_id
        'Montecarlo samples'=mse_df.idx.nunique(), 
        'Bad Samples'=len(bad_ids), 
        'Montecarlo cleaned'=mse_df1.idx.nunique(), 
    }) #.to_csv(os.path.join(out_dpath, f'cleanup_stats_{session_id}.csv')
    

    sim_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_df.csv.gz',), index=False)
    mse_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_mse.csv.gz',), index=False)
    sum_df1.to_csv(os.path.join(out_dpath,f'{session_id}_clean_sum.csv.gz',), index=False)

    return(sim_df1, mse_df1, sum_df1, stats_df)
    
# represent FL 0.1
NBIOMASS_LOD_NOT_GROWING = 2.34502821


def print_vpros(session_id, mse_df, sum_df, stats_df, vprodpath):

    mse_df['RMSE'] = np.sqrt(mse_df['RMSE'])

    mean_scores = mse_df.pivot_table(index=['model', 'idx'], columns='media', values='RMSE' )
    mean_scores.rename(columns={'lowN': 'RMSE_lowN', 'pro99': 'RMSE_pro99'}, inplace=True)
    mean_scores = mean_scores.reset_index()
    mean_scores['RMSE'] = np.sqrt(mean_scores['RMSE_lowN']*  mean_scores['RMSE_pro99'])
    bestids = mean_scores.loc[mean_scores.RMSE.le(60), 'idx'] 


    sum_df1 = sum_df.loc[sum_df.idx.isin(bestids.idx)].drop_duplicates(subset='idx')
    def _get_params_df(model):
        id_vars= ['model', 'idx']
        params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)
        param_vals_df = sum_df1.loc[sum_df1.model.isin([model]), 
                                    id_vars + params_to_update ]
        mparam_vals = param_vals_df.melt( id_vars=id_vars)    
        return mparam_vals
    modellist = sum_df1.model.unique()
    mparams_df = pd.concat([_get_params_df(model) for model in modellist], ignore_index=True)



    for m in modellist:
        os.makedirs(os.path.join(vprodpath, m), exist_ok=True)

    stats_df['VPROs'] = bestids.shape[0]

    for idx in bestids.idx:
        ser_x = mparams_df.loc[
            mparams_df.idx.isin([idx]), ['variable', 'value']
        ]
        ser_x.index= ser_x.variable
        actual_finalX = ser_x.value.to_dict()
        vpro_id = idx.replace(vproprefix,'vpro_').replace('_monte_','')
        fname = vpro_id
        model = mparams_df.loc[mparams_df.idx.isin([idx])]['model'].unique()[0]
        res_fpath = os.path.join(vprodpath, model, f'{fname}.json')
        print(res_fpath)
        params2json(actual_finalX, res_fpath)

    return stats_df




if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='cleanup PONLY runs and produce VPROS.')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--vprooutdpath", help="output dir", default='.')
    parser.add_argument("--indpath", help="input dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    
    
    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    vprodpath = args.vprooutdpath
    if vprodpath != '':
        os.makedirs(vprodpath, exist_ok=True)
    

    dpath = args.indpath
    
    session_id = args.run_id
    sim_df = pd.read_csv(os.path.join(dpath,f'{session_id}_df.csv.gz',))
    mse_df = pd.read_csv(os.path.join(dpath,f'{session_id}_mse.csv.gz',))
    sum_df = pd.read_csv(os.path.join(dpath,f'{session_id}_sum.csv.gz',))

    (sim_df, mse_df, sum_df, stats_df) = cleanup_session(session_id, sim_df, mse_df, sum_df, out_dpath)
    
    stats_df = print_vpros(session_id, mse_df, sum_df, stats_df, vprodpath)
    stats_df.to_csv(os.path.join(out_dpath, f'cleanup_stats_{session_id}.csv')


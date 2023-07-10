#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pprint
sns.set(style="white", context='poster')
from sympy import *
import math
init_printing(use_unicode=True)
from scipy.integrate import solve_ivp
import glob
import os

from model_equations_separate_NC_sep_vmax import *

def _read_csv_try(f):
    try:
        df = pd.read_csv(f, #compression=None
                )    
        return df
    except BaseException as err:
        print(f"Unexpected {err}, {type(err)}")
        print('failed to read:', f)
        return None


def _read_csv_df(f):
    df = _read_csv_try(f)
    if df is not None:
        df['run_id'] = os.path.basename(f) .replace('_df.csv.gz', '') .replace('_mse.csv.gz', '') .replace('_sum.csv.gz', '')
    return df

def create_mse_df(df_fpath, df, refdf):
    if (df is None) or (refdf is None):
        return(0)
    mse_fpath = df_fpath.replace('_df.csv.gz', '_mse.csv.gz')
    if os.path.exists(mse_fpath):
        return(0)
    try:
        mse_df = compute_mse(df, refdf, refcol= 'ref_Bp', col='Bp', timecol='t', tolerance=100)
        mse_df.to_csv(mse_fpath, compression='gzip')
    except BaseException as err:
        print(f"Unexpected {err}, {type(err)}")
        print('failed to create mse:', df_fpath)
        return None


def concat_csvs(dpaths, out_dpath, out_fprefix, refdf):
    res_glob_pattern = '*_df.csv.gz'
    sum_glob_pattern = '*_sum.csv.gz'
    mse_glob_pattern = '*_mse.csv.gz'
    lsq_glob_pattern = '*_least_square.csv'

    print(dpaths)

    data_dfs = [ (f,_read_csv_df(f)) for dpath in dpaths for f in glob.glob(os.path.join(dpath,res_glob_pattern ))] 
    # create mse file if not found
    if refdf is not None:
        [ create_mse_df(f, d, refdf) for f,d in data_dfs if d is not None]
    data_df = pd.concat ( [ d for f,d in data_dfs if d is not None])
    data_df.drop(columns=['Unnamed: 0',], inplace=True)
    data_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'), compression='gzip', index=False)

    sum_dfs = [_read_csv_try(f) for dpath in dpaths for f in glob.glob(os.path.join(dpath,sum_glob_pattern )) ]
    sum_df = pd.concat ( [ d for d in sum_dfs if d is not None])
    sum_df.drop(columns=['Unnamed: 0',], inplace=True)
    sum_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'), compression='gzip', index=False)


    mse_dfs = [ _read_csv_df(f) for dpath in dpaths for f in glob.glob(os.path.join(dpath,mse_glob_pattern ))] 
    if mse_dfs:
        mse_df = pd.concat ( [ d for d in mse_dfs if d is not None])
        mse_df.drop(columns=['Unnamed: 0',], inplace=True)
        mse_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_mse.csv.gz'), compression='gzip', index=False)

    lsq_dfs = [ _read_csv_df(f) for dpath in dpaths for f in glob.glob(os.path.join(dpath,lsq_glob_pattern ))] 
    if lsq_dfs:
        lsq_df = pd.concat ( [ d for d in lsq_dfs if d is not None])
        #lsq_df.drop(columns=['Unnamed: 0',], inplace=True)
        lsq_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_lsq.csv.gz'), compression='gzip', index=False)




def run_umap_hdbscan(df, sum_df):
    import hdbscan
    import umap
    sample_cols = ['Bp', 'Bh', 'DON', 'DIN', 'RDON', 'DOC', 'DIC', 'RDOC']
    sample_days = [ 5., 12., 16., 30., 44., 62.]
    mdf = df.melt(id_vars=['run_id', 'day', 't'], value_vars=ccxorder+ccnorder+cccorder)
    tdf = mdf.loc[mdf.day.round(0).isin(sample_days) & mdf.variable.isin(sample_cols)]
    tdf['day'] = tdf.day.round(0)
    X = tdf.pivot(index='run_id', columns=['variable', 'day'], values=['value'])
    reducer = umap.UMAP(random_state=42)
    umap_embedding = reducer.fit_transform(X)
    udf = pd.DataFrame(data=umap_embedding, columns=['UMAP1','UMAP2'])
    udf['run_id'] = X.index
    udf = pd.merge(udf, sum_df, on='run_id', how='left')

    clusterer = hdbscan.HDBSCAN(min_cluster_size=100, min_samples=15)
    clusterer.fit(X)
    udf['cluster'] = clusterer.labels_
    return udf


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='post process results')
    parser.add_argument('--dpath', #action='append', 
            help='paths to load from', required=True, nargs="+")
    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument('--ref_csv', help='reference CSV', default='None')

    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
        
    if args.ref_csv == 'None':
        refdf = None
    else:
        refdf = pd.read_excel(args.ref_csv)

    concat_csvs(args.dpath, args.outdpath, args.run_id, refdf)

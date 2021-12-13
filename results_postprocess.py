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



def _read_csv_df(f):
    try:
        df = pd.read_csv(f)    
        df['run_id'] = os.path.basename(f).replace('_df.csv.gz', '')
        return df
    except:
        print(f)
        return None

def concat_csvs(dpaths, out_dpath, out_fprefix):
    res_glob_pattern = '*_df.csv.gz'
    sum_glob_pattern = '*_sum.csv.gz'

    sum_dfs = [pd.read_csv(f) for dpath in dpaths for f in glob.glob(os.path.join(dpath,sum_glob_pattern )) ]
    sum_df = pd.concat ( [ d for d in sum_dfs if d is not None])
    sum_df['error'] = sum_df.h_err + sum_df.p_err
    sum_df['logerror'] = np.log(sum_df['error'])
    sum_df.drop(columns=['Unnamed: 0',], inplace=True)
    sum_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'))


    data_dfs = [ _read_csv_df(f) for dpath in dpaths for f in glob.glob(os.path.join(dpath,res_glob_pattern ))] 
    data_df = pd.concat ( [ d for d in data_dfs if d is not None])
    data_df.drop(columns=['Unnamed: 0',], inplace=True)
    data_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'))


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
    parser.add_argument('--dpath', action='append', help='paths to load from', required=True)
    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)

    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
        
    concat_csvs(args.dpath, args.outdpath, args.run_id)

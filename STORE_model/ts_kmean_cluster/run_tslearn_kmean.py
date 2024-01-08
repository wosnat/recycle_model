#!/usr/bin/env python
# coding: utf-8

import itertools
import sys, os

import numpy as np
import pandas as pd
from scipy.special import comb
from scipy import stats
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from tslearn.utils import to_time_series, to_time_series_dataset
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import TimeSeriesKMeans
from tslearn.utils import to_sklearn_dataset

if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Do clusteing.')
    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--model", help="model to run", 
                        choices=['MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME', 'all'], required=True)
    parser.add_argument('--fpath', #action='append', 
            help='paths to load from', required=True, nargs="+")
    parser.add_argument("--classes", help="predicted classes", 
        default='ML_model/monte_predicted_classes.csv.gz')
    parser.add_argument("--n_clusters", help="number of clusters",  type=int, default=3)
    parser.add_argument("--n_init", help="n_init",  type=int, default=1)
    parser.add_argument("--max_iter", help="max_iter",  type=int, default=100)
    parser.add_argument("--max_iter_barycenter", help="max_iter_barycenter",  type=int, default=100)
    parser.add_argument("--n_jobs", help="n_jobs",  type=int, default=-1)
    parser.add_argument("--metric", help="metric",  default='dtw')


    
    
    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)


    model = args.model

    sim_df = pd.concat([pd.read_csv(i) for i in args.fpath], ignore_index=True)

    df_predicted_classes = pd.read_csv(args.classes)

    runids = df_predicted_classes['run_id']
    if model != 'all':
        runids = df_predicted_classes.loc[df_predicted_classes.model.isin([model]), 'run_id']

    sim_df = sim_df.loc[sim_df.run_id.isin(runids), ['run_id', 'day', 'Bptotal[N]'] ].copy()

    
    
    #todo
    X = sim_df.pivot(columns='day', index='run_id', values='Bptotal[N]')
    X_logged = np.log10(X.clip(lower=1e-10))
    tslearn_dataset = to_time_series_dataset(X_logged)


    #tslearn_dataset_norm = TimeSeriesScalerMeanVariance().fit_transform(tslearn_dataset)
    tslearn_dataset_norm = tslearn_dataset

    dba_km = TimeSeriesKMeans(
        n_clusters=args.n_clusters,
        n_init=args.n_init, 
        max_iter=args.max_iter,
        metric=args.metric,
        verbose=True,
        max_iter_barycenter=args.max_iter_barycenter,
        n_jobs=args.n_jobs,
        random_state=1234
        )
    km_y_pred = dba_km.fit_predict(tslearn_dataset_norm)


    cluster_df = pd.DataFrame({'cluster': km_y_pred}, index=X_logged.index)
    cluster_df = cluster_df.reset_index()
    cluster_df = pd.merge(cluster_df, df_predicted_classes, on='run_id', how='left')


    tt = to_sklearn_dataset(tslearn_dataset_norm)
    sample_df = pd.DataFrame(tt, columns=X_logged.columns, index=X_logged.index)
    msample_df = pd.melt(sample_df, ignore_index=False).reset_index()
    msample_df = pd.merge(msample_df, cluster_df, on='run_id')

    df_cluster_centers = pd.DataFrame(to_sklearn_dataset(dba_km.cluster_centers_), 
                                      columns=X_logged.columns)
    mdf_cluster_centers = pd.melt(df_cluster_centers, ignore_index=False).reset_index(names='cluster')

    fprefix = os.path.join(dpath, f'tskmean_{args.run_id}_{model}_{args.n_clusters}')
    dba_km.to_json(f'{fprefix}_kmean.json')
    mdf_cluster_centers.to_csv(f'{fprefix}_cluster_centers.csv.gz')
    msample_df.to_csv(f'{fprefix}_X.csv.gz')
    cluster_df.to_csv(f'{fprefix}_cluster.csv.gz')





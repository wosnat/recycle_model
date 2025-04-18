#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import math
from scipy.stats import gmean
import os
from model_equations_separate_NC_store_numba import *


# adapted from https://github.com/scipy/scipy/blob/v1.11.3/scipy/optimize/_lsq/least_squares.py#L241-L963
# Loss functions.

refdf = pd.read_excel('reference_10cc.xlsx')
refdf = refdf.loc[~refdf.Group.isin(['Axenic'])]

def huber(z):
    rho = np.empty_like(z)
    mask = z <= 1
    rho[mask] = z[mask]
    rho[~mask] = 2 * z[~mask]**0.5 - 1
    return rho

def softl1(z):
    t = 1 + z
    rho = 2 * (t**0.5 - 1)
    return rho

def linear(z):
    t = 1 + z
    rho = 2 * (t**0.5 - 1)
    return rho


def cauchy(z):
    rho = np.log1p(z)
    return rho

def arctan(z):
    rho = np.arctan(z)
    return rho


IMPLEMENTED_LOSSES = dict(linear=linear, huber=huber, softl1=softl1,
                          cauchy=cauchy, arctan=arctan)



def loss_function(f, loss, f_scale):
    lossfunc = IMPLEMENTED_LOSSES[loss]
    z = (f / f_scale) ** 2
    rho = lossfunc(z)
    return 0.5 * f_scale ** 2 * np.sum(rho)

def all_loss_functions(f, f_scale):
    z = (f / f_scale) ** 2
    coef = 0.5 * f_scale ** 2
    return pd.Series({
        f'{loss}_{f_scale}': coef * np.sum(lossfunc(z)) 
        for loss, lossfunc in IMPLEMENTED_LOSSES.items()
    })
    

def _mse_loss(x, refdf, refcol, col, timecol, tolerance, loss, f_scale):
    #print(x.columns)
    tdf = pd.merge_asof(x, refdf[[timecol] + refcol], on=timecol, direction='nearest', tolerance=tolerance).dropna()
    colname = f'{loss}_{f_scale}'
    #mse_series = pd.concat([all_loss_functions(f'RMSE_{c}', tdf[c] - tdf[rc], f_scale) for c,rc in zip(col, refcol)])
    mse_dict = {f'{colname}_{c}': loss_function(tdf[c] - tdf[rc], loss, f_scale) for c,rc in zip(col, refcol)}
    mse_dict[colname] = gmean([mse_dict[f'{colname}_{c}'] for c in col])
    mse_dict['compare_points'] = tdf.shape[0]
    return pd.Series(mse_dict)
    

def compute_mse_loss(
    df, refdf, 
    refcol=['ref_Bp[N]' , 'ref_Bp[C]'], 
    col=['Bptotal[N]', 'Bptotal[C]'], 
    timecol='t', tolerance=100, loss='linear', f_scale=1,
):
    df1 = df[[timecol] + col].dropna().copy()
    df1[col] = df1[col].clip(lower=4)
    mse_df = refdf.groupby(['Sample', 'full name', 'Group',]
                    ).apply(lambda y : 
                            _mse_loss(df1, y, refcol= refcol, col=col, timecol=timecol, tolerance=tolerance, f_scale=f_scale, loss=loss))
    
    mse_df = mse_df.reset_index()
    return mse_df


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='post process results')
    parser.add_argument('--incsv', help='file with results', required=True)
    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument('--ref_csv', help='reference CSV', default='None')
    parser.add_argument('--loss', help='loss function', choices=['linear', 'huber', 'softl1', 'cauchy', 'arctan'], required=True)
    parser.add_argument("--f_scale", help="f_scale",  type=float, required=True)
    parser.add_argument("--which_organism", help="which organism to run", choices=['ponly', 'honly', 'all'], default='all')

    args = parser.parse_args()
    outdpath = args.outdpath
    if outdpath != '':
        os.makedirs(outdpath, exist_ok=True)
        
    if args.ref_csv == 'None':
        refdf = None
    else:
        refdf = pd.read_excel(args.ref_csv)
        refdf['ref_Bp[N]'] = refdf['ref_Bp[N]'].clip(lower=4)
        refdf['ref_Bp[C]'] = refdf['ref_Bp[C]'].clip(lower=4)



    loss = args.loss
    f_scale=args.f_scale
    pro99_mode = False 
    which_organism = args.which_organism
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
    ) = get_constants_per_organism(pro99_mode, which_organism)


    df = pd.read_csv(args.incsv)


    min_values = df[var_names+ ['run_id']].groupby('run_id').min().min(axis=1)
    bad_ids_negative_values = set(min_values[min_values< -1e-9].index)
    bad_ids =bad_ids_negative_values
    df     =     df.loc[    ~df.run_id.isin(bad_ids)].copy()

    msedf = df.groupby(['run_id']).apply(
            lambda x: compute_mse_loss(x, refdf=refdf, f_scale=f_scale, loss=loss)).reset_index()

    msedf.to_csv(os.path.join(outdpath, f'{args.run_id}_loss_{loss}_{f_scale}.csv.gz'), compression='gzip', index=False)


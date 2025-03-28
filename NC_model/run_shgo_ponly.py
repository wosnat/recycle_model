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
from scipy.optimize import least_squares
from functools import lru_cache

from model_equations_separate_NC_sep_vmax import *

from scipy.optimize import minimize
from scipy.optimize import shgo


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', required=True)
    parser.add_argument('--ref_pro99_csv', help='reference pro99 CSV', required=True)
    parser.add_argument('--json_dpath', help='folder to put json files with param vals', default='.')
    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--timeout", help="timeout",  type=int, default=2*60)
    parser.add_argument("--model", help="model to run", choices=['MIN', 'FULL', 'LEAK', 'MIXO'], default='MIN')
    parser.add_argument("--f_scale", help="f_scale",  type=float, default=1)
                        
    
    args = parser.parse_args()
    dpath = args.out_dpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    os.makedirs(args.json_dpath, exist_ok=True)



model = args.model
organism_to_tune = 'PRO'
which_organism='ponly'
run_id = args.run_id

json_dpath = args.json_dpath
out_dpath = args.out_dpath
timeout= args.timeout


ref_fpath =  'reference_10cc_axenic.xlsx'
ref_pro99_fpath =  'refence_pro99_axenic.xlsx'
ref_fpath =  args.ref_csv
ref_pro99_fpath =  args.ref_pro99_csv

ref_df = pd.read_excel(ref_fpath)
ref_pro99_df = pd.read_excel(ref_pro99_fpath)
ref_df = ref_df.sort_values(['t','Sample'])
ref_pro99_df = ref_pro99_df.sort_values(['t','Sample'])
#ref_df = ref_df.loc[ref_df['day'] < 60]
#ref_pro99_df = ref_pro99_df.loc[ref_pro99_df['day'] < 60]
Y = pd.concat([ref_df['ref_Bp'], ref_pro99_df['ref_Bp']]).values
Y = Y.clip(min=4)


param_vals = get_param_vals(model)
params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)

# start with the defalt params
x0 = np.array([np.log(param_vals[i]) if lg else param_vals[i] for i, lg in zip(params_to_update, log_params)])

bounds_logged = [(np.log(b[0]),  np.log(b[1]))  if lg else b for b,lg in zip(bounds, log_params)]
param_bounds =  bounds_logged

t_eval = np.rint(ref_df['t'].drop_duplicates().sort_values()).values
t_eval_pro99 = np.rint(ref_pro99_df['t'].drop_duplicates().sort_values()).values

def run_model(X):
    lowN_run = generate_json_and_run_from_X(
        X, params_to_update, param_vals, 
        'None', json_dpath, out_dpath, run_id, 
            timeout=timeout, log_params=log_params,
        which_organism=which_organism, pro99_mode=False, 
        t_eval = t_eval
    )
    fpath = os.path.join(out_dpath, f'{lowN_run}_df.csv.gz')
    lowN_df = pd.read_csv(fpath)
    result_lowN = pd.merge_asof(ref_df, lowN_df, on='t', tolerance=1, direction='nearest')['Bp']

    pro99_run = generate_json_and_run_from_X(
        X, params_to_update, param_vals, 
        'None', json_dpath, out_dpath, run_id, 
            timeout=timeout, log_params=log_params,
        which_organism=which_organism, pro99_mode=True, 
        t_eval = t_eval_pro99
    )
    fpath = os.path.join(out_dpath, f'{pro99_run}_df.csv.gz')
    pro99_df = pd.read_csv(fpath)
    result_pro99 = pd.merge_asof(ref_pro99_df, pro99_df, on='t', tolerance=1, direction='nearest')['Bp']
    return pd.concat([result_lowN, result_pro99]).values.clip(min=4)

# this wrap is to cache run_model results
run_model = lru_cache(run_model)
def wrap_run_model(X):
    return(run_model(tuple(X)))




def jac(X):
    X = np.array(X)
    #print(X)
    delta = 1e-8
    base = wrap_run_model(X)
    J = np.empty((base.size, X.size))
    for i in range(X.size):
        deltaX = np.array(X)
        deltaX[i] = deltaX[i] +  delta
        delta_Bp = wrap_run_model(deltaX)
        J[:, i] = (delta_Bp - base) / delta

    return J    

jac = lru_cache(jac)
def wrap_jac(X):
    return(jac(tuple(X)))


def fun(X):
    err = wrap_run_model(X) - Y
    return np.sum(err**2)



def jac_for_minimize(X):
    J = wrap_jac(X)
    err = wrap_run_model(X) - Y
    return J.T.dot(err)


def hess_for_minimize(X):
    J = wrap_jac(X)
    return J.T.dot(J)




print('shgo(',
    fun, param_bounds, 
    #minimizer_kwargs=dict(
        #jac=jac_for_minimize, 
#        hess=hess_for_minimize, 
#        options=dict(disp=True, maxiter=3), method='L-BFGS-B'),
#    options=dict(jac=
            jac_for_minimize)

res = shgo(
    fun, param_bounds, 
    n=100, iters=1, sampling_method='sobol',
    minimizer_kwargs=dict(
        #jac=jac_for_minimize, 
        hess=hess_for_minimize, 
        options=dict(disp=True, maxiter=3), method='L-BFGS-B'),
    options=dict(jac=jac_for_minimize,disp=True),
)

print(res)

for idx, finalX in enumerate(res.xl, start=1):
    actual_finalX = {p: np.exp(i) if lg else i for i,lg,p in zip(finalX, log_params, params_to_update)}
    res_fpath = os.path.join(out_dpath, f'{run_id}_{idx}.json')
    params2json(actual_finalX, res_fpath)





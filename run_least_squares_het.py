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


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', required=True)
    parser.add_argument('--ref_sample_id', help='reference sample id for fitting', required=True, type=int)
    parser.add_argument('--vpro_json', help='json file used for pro params', required=True)
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
organism_to_tune = 'HET'
which_organism='all'
run_id = args.run_id

json_dpath = args.json_dpath
out_dpath = args.out_dpath
timeout= args.timeout
vpro_json = args.vpro_json

ref_fpath =  args.ref_csv
ref_sample_id = args.ref_sample_id

ref_df = pd.read_excel(ref_fpath)
ref_df = ref_df.loc[ref_df['id'].isin([ref_sample_id])]
ref_df = ref_df.sort_values(['t','Sample'])
Y = ref_df['ref_Bp'].values
Y = Y.clip(min=4)


param_vals = get_param_vals(model)
param_vals = json2params(param_vals, vpro_json)

params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)

# start with the defalt params
x0 = np.array([np.log(param_vals[i]) if lg else param_vals[i] for i, lg in zip(params_to_update, log_params)])

bounds_logged = [(np.log(b[0]),  np.log(b[1]))  if lg else b for b,lg in zip(bounds, log_params)]
param_bounds =  list(zip(*bounds_logged))

t_eval = np.rint(ref_df['t'].drop_duplicates().sort_values()).values


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

    return result_lowN.values.clip(min=4)

# this wrap is to cache run_model results
run_model = lru_cache(run_model)
def wrap_run_model(X):
    return(run_model(tuple(X)))

#def jac(X):
#    delta = 1e-8
#    base = wrap_run_model(X)
#    J = np.empty((base.size, X.size))
#    for i in range(X.size):
#        delta_i = X[i]*delta
#        deltaX = np.array(X)
#        deltaX[i] = deltaX[i] +  delta_i
#        delta_Bp = wrap_run_model(deltaX)
#        J[:, i] = (delta_Bp - base) / delta_i
#
#    return J    
    
def jac(X):
    delta = 1e-8
    base = wrap_run_model(X)
    J = np.empty((base.size, X.size))
    for i in range(X.size):
        deltaX = np.array(X)
        deltaX[i] = deltaX[i] +  delta
        delta_Bp = wrap_run_model(deltaX)
        J[:, i] = (delta_Bp - base) / delta

    return J    

def fun(X):
    return wrap_run_model(X) - Y

f_scale=args.f_scale
res = least_squares(fun, x0, jac=jac, bounds=param_bounds,  verbose=2,  loss='soft_l1', f_scale=f_scale)
print(res)

finalX = res.x
actual_finalX = {p: np.exp(i) if lg else i for i,lg,p in zip(finalX, log_params, params_to_update)}


res_fpath = os.path.join(out_dpath, f'{run_id}.json')
params2json(actual_finalX, res_fpath)





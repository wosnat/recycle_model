#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pprint
sns.set(style="white", context='poster')
import math
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from functools import lru_cache

from model_equations_separate_NC_store_numba import *


def run_model(X, additional_params):
    (Y, params_to_update, orig_param_vals, log_params, 
    calc_dydt, prepare_params_tuple, var_names, 
    init_var_vals, t_end , t_eval, ref_df, 
    init_var_pro99_vals, t_end_pro99 , t_eval_pro99, ref_pro99_df, logerror
    ) = additional_params
    #print(X)

    try:
        new_param_vals = get_params(X, params_to_update, orig_param_vals, log_params)
        par_tuple = prepare_params_tuple(new_param_vals)

        lowN_sol = run_solver(calc_dydt, init_var_vals, par_tuple, t_end , t_eval)
        lowN_df = solver2df_forlsq(lowN_sol, var_names, par_tuple)
        result_lowN = pd.merge_asof(ref_df, lowN_df, on='t', tolerance=1, direction='nearest')

        pro99_sol = run_solver(calc_dydt, init_var_pro99_vals, par_tuple, t_end_pro99 , t_eval_pro99)
        pro99_df = solver2df_forlsq(pro99_sol, var_names, par_tuple)
        result_pro99 = pd.merge_asof(ref_pro99_df, pro99_df, on='t', tolerance=1, direction='nearest')

        res = np.clip(pd.concat([
                result_lowN['Bptotal[N]'], 
                result_lowN['Bptotal[C]'], 
                result_pro99['Bptotal[N]'],
                result_pro99['Bptotal[C]'],
            ]).to_numpy(), a_min=4, a_max=None)
        if logerror:
            res = np.log(res)
        return res
    except BaseException as err:
        print(f"Unexpected {err}, {type(err)}")
        return np.zeros_like(Y)

def jac(X, additional_params):
    delta = 1e-8
    base = run_model(X, additional_params)
    J = np.empty((base.size, X.size))
    for i in range(X.size):
        deltaX = np.array(X)
        deltaX[i] = deltaX[i] +  delta
        delta_Bp = run_model(deltaX, additional_params)
        J[:, i] = (delta_Bp - base) / delta

    return J    

def fun(X, additional_params):
    return run_model(X, additional_params) - Y


if __name__ == '__main__':
    import argparse
    import json
    import pprint
    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', required=True)
    parser.add_argument('--ref_pro99_csv', help='reference pro99 CSV', required=True)
    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--model", help="model to run", choices=['MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME'], required=True)
    parser.add_argument('--json', help='json with param vals', nargs="+")
    parser.add_argument('--logerror', help='use log of error', action="store_true")
                        
    
    args = parser.parse_args()
    dpath = args.out_dpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)

    model = args.model
    organism_to_tune = 'PRO'
    which_organism='ponly'
    run_id = args.run_id

    out_dpath = args.out_dpath


    ref_fpath =  args.ref_csv
    ref_pro99_fpath =  args.ref_pro99_csv

    ref_df = pd.read_excel(ref_fpath)
    ref_pro99_df = pd.read_excel(ref_pro99_fpath)
    ref_df = ref_df.sort_values(['t','Sample'])
    ref_pro99_df = ref_pro99_df.sort_values(['t','Sample'])

    Y = pd.concat([
            ref_df['ref_Bp[N]'], 
            ref_df['ref_Bp[C]'], 
            ref_pro99_df['ref_Bp[N]']
            ref_pro99_df['ref_Bp[C]']
        ]).to_numpy()
    Y = np.clip(Y, a_min=4, a_max=None)
    if args.logerror:
        Y = np.log(Y)

    new_param_vals = get_param_vals_from_json_list(args.model, args.json)
    #TODO
    t_eval, t_end = get_t_eval_and_t_end(None, ref_df, maxday=140)
    t_eval_pro99, t_end_pro99 = get_t_eval_and_t_end(None, ref_pro99_df, maxday=140)
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
        ) = get_constants_per_organism(False, which_organism)
    (_, init_var_pro99_vals, _, _, _) = get_constants_per_organism(True, which_organism)


    params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)

    # start with the defalt params
    x0 = np.array([np.log(new_param_vals[i]) if lg else new_ [i] for i, lg in zip(params_to_update, log_params)])

    bounds_logged = [(np.log(b[0]),  np.log(b[1]))  if lg else b for b,lg in zip(bounds, log_params)]
    param_bounds =  list(zip(*bounds_logged))

    additional_params = (
        Y, params_to_update, new_param_vals, log_params, 
        calc_dydt, prepare_params_tuple, var_names, 
        init_var_vals, t_end , t_eval, ref_df, 
        init_var_pro99_vals, t_end_pro99 , t_eval_pro99, ref_pro99_df, args.logerror
    ) 

    f_scale=1
    res = least_squares(fun, x0, jac=jac, bounds=param_bounds,  verbose=2,  loss='soft_l1', f_scale=f_scale, args=(additional_params,))
    print(res)

    finalX = res.x
    actual_finalX = {p: np.exp(i) if lg else i for i,lg,p in zip(finalX, log_params, params_to_update)}


    res_fpath = os.path.join(out_dpath, f'{run_id}.json')
    params2json(actual_finalX, res_fpath)





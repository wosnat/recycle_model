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
from scipy.stats import qmc
from functools import lru_cache

from model_equations_separate_NC_store_numba import *

from scipy.optimize import minimize
from scipy.optimize import shgo






def run_model(X, additional_params):
    (Y, params_to_update, orig_param_vals, log_params, 
    calc_dydt, prepare_params_tuple, var_names, 
    init_var_vals, t_end , t_eval, ref_df, 
    init_var_pro99_vals, t_end_pro99 , t_eval_pro99, ref_pro99_df
    ) = additional_params
    #print(X)

    try:
        new_param_vals = get_params(X, params_to_update, orig_param_vals, log_params)
        par_tuple = prepare_params_tuple(new_param_vals)

        lowN_sol = run_solver(calc_dydt, init_var_vals, par_tuple, t_end , t_eval)
        lowN_df = solver2df_forlsq(lowN_sol, var_names)
        result_lowN = pd.merge_asof(ref_df, lowN_df, on='t', tolerance=1, direction='nearest')['Bptotal']

        pro99_sol = run_solver(calc_dydt, init_var_pro99_vals, par_tuple, t_end_pro99 , t_eval_pro99)
        pro99_df = solver2df_forlsq(pro99_sol, var_names)
        result_pro99 = pd.merge_asof(ref_pro99_df, pro99_df, on='t', tolerance=1, direction='nearest')['Bptotal']

        return np.clip(pd.concat([result_lowN, result_pro99]).to_numpy(), a_min=4, a_max=None)
    except BaseException as err:
        print(f"Unexpected {err}, {type(err)}")
        return np.zeros_like(Y)

# this wrap is to cache run_model results
# run_model = lru_cache(run_model)
#def wrap_run_model(X, additional_params):
#    return(run_model(tuple(X)))

def jac(X, additional_params):
    X = np.array(X)
    #print(X)
    delta = 1e-8
    base = run_model(X, additional_params)
    J = np.empty((base.size, X.size))
    for i in range(X.size):
        deltaX = np.array(X)
        deltaX[i] = deltaX[i] +  delta
        delta_Bp = run_model(deltaX, additional_params)
        J[:, i] = (delta_Bp - base) / delta

    return J    

# jac = lru_cache(jac)
# def wrap_jac(X):
    # return(jac(tuple(X)))


def fun(X, additional_params):
    Y = additional_params[0]
    err = run_model(X, additional_params) - Y
    return np.sum(err**2)



def jac_for_minimize(X, additional_params):
    Y = additional_params[0]
    J = jac(X, additional_params)
    err = run_model(X, additional_params) - Y
    return J.T.dot(err)


def hess_for_minimize(X, additional_params):
    J = jac(X, additional_params)
    return J.T.dot(J)

def get_sobol_sample(param_bounds, m):
    ''' return a list of samples
    param_bounds - list of (lower,upper) bounds for scaling
    m - will produce 2^m samples
    '''
    l_bounds, u_bounds = tuple(zip(*param_bounds))
    sampler = qmc.Sobol(d=len(param_bounds))
    sample = sampler.random_base2(m=m)
    return qmc.scale(sample, l_bounds, u_bounds)




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
    parser.add_argument("--number_of_runs", help="number of simulations to run",  type=int, default=1024)
                        
    
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

    Y = pd.concat([ref_df['ref_Bp'], ref_pro99_df['ref_Bp']]).to_numpy()
    Y = np.clip(Y, a_min=4, a_max=None)

    new_param_vals = get_param_vals_from_json_list(args.model, args.json)
    #TODO
    t_eval, t_end = get_t_eval_and_t_end(None, ref_df, maxday=140)
    t_eval_pro99, t_end_pro99 = get_t_eval_and_t_end(None, ref_pro99_df, maxday=140)
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
        ) = get_constants_per_organism(False, which_organism)
    (_, init_var_pro99_vals, _, _, _) = get_constants_per_organism(True, which_organism)


    params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)

    # start with the defalt params
    x0 = np.array([np.log(param_vals[i]) if lg else param_vals[i] for i, lg in zip(params_to_update, log_params)])

    bounds_logged = [(np.log(b[0]),  np.log(b[1]))  if lg else b for b,lg in zip(bounds, log_params)]
    param_bounds =  bounds_logged

    # print('shgo(',
        # fun, param_bounds, 
        # #minimizer_kwargs=dict(
            # #jac=jac_for_minimize, 
    # #        hess=hess_for_minimize, 
    # #        options=dict(disp=True, maxiter=3), method='L-BFGS-B'),
    # #    options=dict(jac=
                # jac_for_minimize)

    additional_params = (
        Y, params_to_update, new_param_vals, log_params, 
        calc_dydt, prepare_params_tuple, var_names, 
        init_var_vals, t_end , t_eval, ref_df, 
        init_var_pro99_vals, t_end_pro99 , t_eval_pro99, ref_pro99_df
    ) 


    sample = get_sobol_sample(param_bounds, m=args.number_of_runs)
    SSE_list = [fun(X, additional_params) for X in sample]
    sample_unlogged = np.copy(sample)
    sample_unlogged[:, log_params] = np.exp(sample_unlogged[:, log_params])
    df = pd.DataFrame(sample_unlogged, columns=params_to_update)
    df['SSE'] = SSE_list          
    #df.nsmallest(n=3, columns='SSE')
    res_fpath = os.path.join(out_dpath, f'{run_id}_sobol_ponly.csv.gz')
    df.to_csv(res_fpath, compression='gzip', index=False)
    
    # can't get SHGO to work (scipy bugs)
    # res = shgo(
        # fun, param_bounds, args=(additional_params,), 
        # n=args.number_of_runs, iters=1, sampling_method='sobol',
        # minimizer_kwargs=dict(
            # jac=jac_for_minimize, 
            # hess=hess_for_minimize, 
            # options=dict(disp=True, maxiter=3, ), method='L-BFGS-B'),
        # options=dict(jac=jac_for_minimize,disp=True),
    # )

    # print(res)

    # for idx, finalX in enumerate(res.xl, start=1):
        # actual_finalX = {p: np.exp(i) if lg else i for i,lg,p in zip(finalX, log_params, params_to_update)}
        # res_fpath = os.path.join(out_dpath, f'{run_id}_{idx}.json')
        # params2json(actual_finalX, res_fpath)


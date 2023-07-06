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


def run_model(random_params, random_values):
    lowN_run = generate_json_and_run_from_X(
        random_values, random_params, param_vals, 
        ref_fpath, json_dpath, out_dpath, run_id, 
            timeout=timeout, log_params=log_params,
        which_organism=which_organism, pro99_mode=False,
    )

# this wrap is to cache run_model results
run_model = lru_cache(run_model)

def wrap_run_model(X):
    return(run_model(tuple(X)))



def create_random_param_vals(i):
    random_number_of_parameters = i[0]
    random_parameter_values = i[1:]
    random_params_to_update = rng.choice(
        len(params_to_update), replace=False, size=random_number_of_parameters)
    random_params = [params_to_update[j] for j in random_params_to_update ]
    random_values = [random_parameter_values[j] for j in random_params_to_update ]
    return random_params, random_values

    
def run_monte_carlo(number_of_runs):
    rng = np.random.default_rng()
    number_of_params = rng.integers(low=2,high=6, size=number_of_runs)
    random_param_values = [rng.uniform(low=l, high=h, size=number_of_runs) for l,h in bounds_logged]
    for i in zip(number_of_params, *random_param_values):
        random_params, random_values = create_random_param_vals(i)
        print (random_params, random_values)
            


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', required=True)
    parser.add_argument('--vpro_json', help='json file used for pro params', required=True)
    parser.add_argument('--json_dpath', help='folder to put json files with param vals', default='.')
    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--timeout", help="timeout",  type=int, default=2*60)
    parser.add_argument("--model", help="model to run", choices=['MIN', 'FULL', 'LEAK', 'MIXO'], default='MIN')
    parser.add_argument("--number_of_runs", help="number of simulations to run",  type=int, default=1024)
          
    
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
    number_of_runs = args.number_of_runs


    param_vals = get_param_vals(model)
    param_vals = json2params(param_vals, vpro_json)

    params_to_update, bounds, log_params = get_param_tuning_values(model, organism_to_tune)

    # start with the defalt params
    # x0 = np.array([np.log(param_vals[i]) if lg else param_vals[i] for i, lg in zip(params_to_update, log_params)])

    bounds_logged = [(np.log(b[0]),  np.log(b[1]))  if lg else b for b,lg in zip(bounds, log_params)]
    param_bounds =  list(zip(*bounds_logged))

    run_monte_carlo(number_of_runs)
#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pprint
sns.set(style="white", context='poster')
from sympy import *
import math
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error
import json
import subprocess
import sys
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pprint
from sympy import *
import math
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution


from SALib.sample import saltelli
from SALib.analyze import sobol

from model_equations_separate_NC_sep_vmax import *

def generate_samples():

    params_to_update = [
        'M_h', 'M_p', 'gamma^D_p', 'gamma^D_h', 
        #'R_p', 'R_h', 
        'E^O_p', 'E^I_p', 'E^O_h', 'E^I_h', 
        'K^ON_p', 'K^IN_p', 'K^OC_p', 'K^IC_p', 'K^ON_h', 'K^IN_h', 'K^OC_h', 'K^IC_h', 
        # 'Vmax^ON_p', 
        'Vmax^IN_p', 
        # 'Vmax^OC_p', 
        'Vmax^IC_p', 'Vmax^ON_h', 'Vmax^IN_h', 'Vmax^OC_h', 
        #'Vmax^IC_h', 
        'O_p', 'O_h', 'epsilon', 'VTmax', 'KT_h', 'omega'
    ]

    bounds = [
        (0/ seconds_in_day,1/ seconds_in_day ), # 'M_h', 
        (0/ seconds_in_day,1/ seconds_in_day ), # 'M_p', 
        (0,1 ), # 'gamma^D_p', 
        (0,1 ), # 'gamma^D_h', 
        # ?? 'R_p', 'R_h', 
        
        (0/ seconds_in_day,0.9/ seconds_in_day ), # 'E^O_p', 
        (0/ seconds_in_day,0.9/ seconds_in_day ), # 'E^I_p',
        (0/ seconds_in_day,0.9/ seconds_in_day ), # 'E^O_h',
        (0/ seconds_in_day,0.9/ seconds_in_day ), # 'E^I_h',
        (param_vals[str(KONp)]/100,param_vals[str(KONp)]*100 ), # 'K^ON_p', 
        (param_vals[str(KINp)]/100,param_vals[str(KINp)]*100 ), # 'K^IN_p', 
        (param_vals[str(KOCp)]/100,param_vals[str(KOCp)]*100 ), # 'K^OC_p', 
        (param_vals[str(KICp)]/100,param_vals[str(KICp)]*100 ), # 'K^IC_p', 
        (param_vals[str(KONh)]/100,param_vals[str(KONh)]*100 ), # 'K^ON_h',  
        (param_vals[str(KINh)]/100,param_vals[str(KINh)]*100 ), # 'K^IN_h', 
        (param_vals[str(KOCh)]/100,param_vals[str(KOCh)]*100 ), # 'K^OC_h', 
        (param_vals[str(KICh)]/100,param_vals[str(KICh)]*100 ), # 'K^IC_h', 
        # 'Vmax^ON_p',
        (param_vals[str(VmaxINp)]/100,param_vals[str(VmaxINp)]*100 ), #  'Vmax^IN_p',
        # 'Vmax^OC_p', 
        (param_vals[str(VmaxICp)]/100,param_vals[str(VmaxICp)]*100 ), # 'Vmax^IC_p', 
        (param_vals[str(VmaxONh)]/100,param_vals[str(VmaxONh)]*100 ), # 'Vmax^ON_h', 
        (param_vals[str(VmaxINh)]/100,param_vals[str(VmaxINh)]*100 ), # 'Vmax^IN_h', 
        (param_vals[str(VmaxOCh)]/100,param_vals[str(VmaxOCh)]*100 ), # 'Vmax^OC_h', 
        # 'Vmax^IC_h', 
        (0, 1 ), # 'O_p', 
        (0, 1 ), # 'O_h', 
        (0/ seconds_in_day,1/ seconds_in_day ), # 'epsilon', 
        (param_vals[str(VTmax)]/100,param_vals[str(VTmax)]*100 ), # 'VTmax', 
        (param_vals[str(KTh)]/100,param_vals[str(KTh)]*100 ), # 'KT_h'
        (0,2), # 'omega'
        
    ]

    problem = {
        'num_vars': len(params_to_update),
        'names': params_to_update,
        'bounds': bounds,
    }




if __name__ == '__main__':
    import argparse
    import json
    import pprint

    # def generate_json_and_run_from_X(X,  json_dpath, out_dpath, out_fprefix, timeout=10*60):
    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', default='prelim bottle.csv')
    parser.add_argument('--params_txt', help='parameters file', default='param_values.txt.gz')
    parser.add_argument('--json_dpath', help='folder to put json files with param vals', default='.')
    # parser.add_argument('--workers', help='number of workers', type=int, default=-1)

    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--timeout", help="timeout",  type=int, default=2*60)
    parser.add_argument("--number_of_runs", help="number of simulations to run",  type=int, default=1000)
    parser.add_argument("--chunk", help="which of the chunks to run ",  type=int, default=1)
    parser.add_argument("--optimize", help="run optimization (default: run model)",
                        action="store_true")
    parser.add_argument("--gen_sensitivity", help="generate sensitivity",
                        action="store_true")
    
    args = parser.parse_args()
    dpath = args.out_dpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    os.makedirs(args.json_dpath, exist_ok=True)
    
    refdf = pd.read_csv(args.ref_csv)
    model_name = args.run_id

    params_to_update = [
        'M_h', 'M_p', 'gamma^D_p', 'gamma^D_h', 
        #'R_p', 'R_h', 
        'E^O_p', 'E^I_p', 'E^O_h', 'E^I_h', 
        'K^ON_p', 'K^IN_p', 'K^OC_p', 'K^IC_p', 'K^ON_h', 'K^IN_h', 'K^OC_h', 'K^IC_h', 
        # 'Vmax^ON_p', 
        'Vmax^IN_p', 
        # 'Vmax^OC_p', 
        'Vmax^IC_p', 'Vmax^ON_h', 'Vmax^IN_h', 'Vmax^OC_h', 
        #'Vmax^IC_h', 
        'O_p', 'O_h', 'epsilon', 'VTmax', 'KT_h', 'omega'
    ]


    if args.gen_sensitivity:
        problem = {
            'num_vars': len(params_to_update),
            'names': params_to_update,
            'bounds': bounds,
        }
        param_values = saltelli.sample(problem, 1024)
        print(param_values.shape)
        np.savetxt(args.params_txt, param_values)

    elif args.optimize:
        func = lambda X :  generate_json_and_run_from_X(
            X, params_to_update, param_vals, 
            args.ref_csv, args.jsondpath, dpath, run_id, 
            timeout=args.timeout)
        result = differential_evolution(func, bounds, disp=True, init='latinhypercube')
        print(result.message)
        print(result.x, result.fun)
        print(zip (params_to_update, result.x))
        res_dict = {
            'fun': result.fun,
             'message': result.message,
             'nfev': result.nfev,
             'nit': result.nit,
             'success': result.success,
            'run_id' : args.run_id,
        }
        res_dict.update({i:v for i,v in zip (params_to_update, result.x)})
        pd.DataFrame([res_dict]).to_csv(os.path.join(dpath, f'{args.run_id}_differential_evolution.csv.gz'))
        
    else:
    param_values = np.loadtxt(args.params_txt)
    run_chunk(
        param_values, params_to_update, args.chunk, args.number_of_runs, 
        args.run_id, args.ref_csv, args.json_dpath, args.out_dpath, args.timeout
    )




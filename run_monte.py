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
from itertools import combinations


from SALib.sample import saltelli
from SALib.analyze import sobol

from model_equations_separate_NC_sep_vmax import *


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    # def generate_json_and_run_from_X(X,  json_dpath, out_dpath, out_fprefix, timeout=10*60):
    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', default='prelim bottle.csv')
    parser.add_argument('--json_dpath', help='folder to put json files with param vals', default='.')

    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--timeout", help="timeout",  type=int, default=2*60)
    parser.add_argument("--number_of_runs", help="number of simulations to run",  type=int, default=1000)
    parser.add_argument("--model", help="model to run", choices=['MIN', 'FULL', 'LEAK', 'MIXO'], default='FULL')
    parser.add_argument("--organism_to_tune", help="which organism to tune", choices=['PRO', 'HET'], default='PRO')
    parser.add_argument("--which_organism", help="which organism to run", choices=['ponly', 'honly', 'all'], default='all')
    parser.add_argument('--ref_sample_id', help='reference sample id for fitting', required=True, type=int)
    parser.add_argument('--vpro_json', help='json file used for pro params', required=True)
                        
    
    args = parser.parse_args()
    dpath = args.out_dpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    os.makedirs(args.json_dpath, exist_ok=True)
    
    #refdf = pd.read_excel(args.ref_csv)
    param_vals = get_param_vals(args.model)

    params_to_update, bounds, log_params = get_param_tuning_values(args.model, args.organism_to_tune)
    run_id = f"{args.run_id}_{args.model}"


    if args.gen_sensitivity:
        problem = {
            'num_vars': len(params_to_update),
            'names': params_to_update,
            'bounds': bounds,
        }
        param_values = saltelli.sample(problem, 1024)
        print(param_values.shape)
        np.savetxt(args.params_txt, param_values)

        
    elif args.param_sensitivity != -1:
        idx = args.param_sensitivity
        run_sensitivity_per_parameter(param_vals, params_to_update[idx], bounds[idx], args.number_of_runs, 
            run_id, args.ref_csv, args.json_dpath, args.out_dpath, args.timeout, log_param=log_params[idx], 
            which_organism=args.which_organism, pro99_mode=args.pro99_mode,
        )
        
    else:
        param_values = np.loadtxt(args.params_txt)
        run_chunk(param_vals, 
            param_values, params_to_update, args.chunk, args.number_of_runs, 
            run_id, args.ref_csv, args.json_dpath, args.out_dpath, 
            args.timeout, log_params=log_params, which_organism=args.which_organism, pro99_mode=args.pro99_mode
        )




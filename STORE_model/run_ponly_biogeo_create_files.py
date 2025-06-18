#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", context='poster')
import math
from scipy import integrate

from model_equations_separate_NC_store_numba import *


def _integrate_column(x, collist):
    return pd.Series({colname :  integrate.simpson(x[colname], x=x['t']) for colname in collist})
        
def integrate_variables(df):
    collist = [
        'Bp', 'Np', 'Cp',  'DON', 'RDON', 'DIN',
           'DOC', 'RDOC', 'DIC', 'ROS', 'gross_uptakeINp',
            'gross_uptakeONp', 
           'gross_uptakeICp',  'gross_uptakeOCp',
        'respirationCp',
           'overflowNp', 'overflowCp', 
           'Bptotal[N]', 'Bptotal[C]', 
           'ROSproductionp', 'ROSlossp', 
           'DON2DIN', 'additionalLossRatep', 
        'deathC_DOCp',  'deathN_DONp',
    ]
    total_gross_uptake = df.groupby(['run_id', ]).apply(lambda x: _integrate_column(x, collist))
    total_gross_uptake = total_gross_uptake.reset_index()


    N_uptake_cols = ['gross_uptakeINp', 'gross_uptakeONp', ]
    total_gross_uptake['PP'] = total_gross_uptake['gross_uptakeICp']
    total_gross_uptake['Total N uptake'] = total_gross_uptake[N_uptake_cols].sum(axis=1)
    return total_gross_uptake


##################################
if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='cleanup PONLY runs and produce VPROS.')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--indpath", help="input dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    
    
    args = parser.parse_args()
    out_dpath = args.outdpath
    if out_dpath != '':
        os.makedirs(out_dpath, exist_ok=True)

    dpath = args.indpath
    
    session_id = args.run_id
    sim_df = pd.read_csv(os.path.join(dpath,f'{session_id}_clean_df.csv.gz',))
    

    df = sim_df

    df['total C biomass'] = df['Bptotal[C]'] 
    df['total fixed C'] = df['total C biomass'] + df['DOC'] + df['RDOC']
    


    # # integrate

    total_gross_uptake = integrate_variables(df)
    last_day_df = df.loc[df.groupby('run_id').day.idxmax()]



    id_keys = ['run_id', 
              ]
    comb_df = pd.merge(
        last_day_df, 
        total_gross_uptake, 
        on=id_keys,
        suffixes=[' Final', ' Integrated'],
    )


    comb_df['PP / total fixed C'] = comb_df['PP'] / comb_df['total fixed C'] 
    comb_df['N reuse'] = comb_df['Total N uptake'] / (INIT_DIN + INIT_DON)

    comb_df.to_csv(os.path.join(out_dpath, f'biogeo_{session_id}.csv.gz'), index=False)


#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
from scipy import integrate

from model_equations_separate_NC_store_numba import *

dpath = 'results/final/het'

fpath = 'core_models_versatile_het_df.csv.gz'
df = pd.read_csv(os.path.join(dpath, fpath)) 



def _integrate_column(x, collist):
    return pd.Series({colname :  integrate.simpson(x[colname], x=x['t']) for colname in collist})
        


collist = [
    'Bp', 'Np', 'Cp', 'Bh', 'Nh', 'Ch', 'DON', 'RDON', 'DIN',
       'DOC', 'RDOC', 'DIC', 'ROS', 'gross_uptakeINp',
       'gross_uptakeINh', 'gross_uptakeONp', 'gross_uptakeONh',
       'gross_uptakeICp', 'gross_uptakeICh', 'gross_uptakeOCp',
       'gross_uptakeOCh', 
        #'biosynthesisNp', 'biosynthesisNh', 
    'respirationCp',
       'respirationCh', #'biomass_breakdownCp', 'biomass_breakdownCh',
       'overflowNp', 'overflowNh', 'overflowCp', 'overflowCh', #'Bp[C]',
       'Bptotal[N]', 'Bptotal[C]', #'Bh[C]', 
    'Bhtotal[N]', 'Bhtotal[C]',
       'ROSproductionp', 'ROSproductionh', 'ROSlossp', 'ROSlossh',
       #'deathbiomassNp', 'deathbiomassNh', 'deathstoreNp', 'deathstoreNh',
       #'deathstoreCp', 'deathstoreCh', 
    #'DON2DIN_exop', 
    'DON2DIN_exoh',
       'DON2DIN', 'additionalLossRatep', 'additionalLossRateh', 
    #'regNp',
      # 'regNh', 'regCp', 'regCh', 
    'deathC_DOCp', 'deathC_DOCh', 'deathN_DONp',
       'deathN_DONh',
]


idvars = [
    'model', 'Phase', 'VPRO', 'run_id', 'y_pred', 
   'max_prob', 
]
   
total_gross_uptake = df.groupby(idvars).apply(lambda x: _integrate_column(x, collist))

total_gross_uptake = total_gross_uptake.reset_index()
        
total_gross_uptake.to_csv(os.path.join(dpath, 'integrated_core_models_versatile_het.csv.gz'), index=False)



# def gen_df(y, sum_df=rsum_df, which_organism='all'):
    # run_id = y.name
    # param_vals = sum_df.loc[sum_df['run_id'].isin([run_id])].squeeze().to_dict()
    # (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
        # ) = get_constants_per_organism(pro99_mode=False, which_organism=which_organism)
    # par_tuple = prepare_params_tuple(param_vals)
    # d = y.copy()
    # d[intermediate_names] = d[var_names].apply(
        # lambda x : calc_dydt(0, x.to_numpy(), par_tuple=par_tuple, return_intermediate=True), axis=1, 
        # result_type='expand')
    # return d
    
# df = sim_df1.groupby('run_id',group_keys=True).apply(gen_df)
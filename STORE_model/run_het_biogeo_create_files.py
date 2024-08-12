#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", context='poster')
import math

from model_equations_separate_NC_store_numba import *

dpath = 'results/final/het'


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='cleanup PONLY runs and produce VPROS.')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--infpath", help="input dir", default='.')
    
    
    args = parser.parse_args()
    out_dpath = args.outdpath
    if out_dpath != '':
        os.makedirs(out_dpath, exist_ok=True)

#######################

    fpath = args.infpath
    # results/final/het/monte_het_add_100per_vpro_ROS_clean_df.csv.gz
    fname = os.path.basename(fpath)
    session_id = fname.replace('_df.csv.gz', '')
    
    df = pd.read_csv(fpath)
    df = df.reset_index(drop=True)



fnames_predicted_classes = [
    'predicted_classes_monte_add_het_clean.csv.gz',
    'predicted_classes_monte_het_clean.csv.gz',
    'predicted_classes_monte_het_add_100per_vpro_EXOENZYME.csv.gz',
    'predicted_classes_monte_het_add_100per_vpro_OVERFLOW.csv.gz',
    'predicted_classes_monte_het_add_100per_vpro_ROS.csv.gz',
    'predicted_classes_monte_het_extend_100per_vpro_EXOENZYME.csv.gz',
    'predicted_classes_monte_het_extend_100per_vpro_MIXOTROPH.csv.gz',
    'predicted_classes_monte_het_extend_100per_vpro_OVERFLOW.csv.gz',
    'predicted_classes_monte_het_extend_100per_vpro_ROS.csv.gz',
    'predicted_classes_monte_het_extend_EXOENZYME.csv.gz',
    'predicted_classes_monte_het_extend_MIXOTROPH.csv.gz',
    'predicted_classes_monte_het_extend_OVERFLOW.csv.gz',
    'predicted_classes_monte_het_extend_ROS.csv.gz',
    'predicted_classes_monte_het_multi.csv.gz',
    'predicted_classes_monte_het_round2_100per_vpro_ROS.csv.gz',
    'predicted_classes_monte_ROS_round2_het.csv.gz',
]


# In[ ]:


fnames_df = [
    'monte_add_het_clean_clean_df.csv.gz',
    'monte_het_add_100per_vpro_EXOENZYME_clean_df.csv.gz',
    'monte_het_add_100per_vpro_OVERFLOW_clean_df.csv.gz',
    'monte_het_add_100per_vpro_ROS_clean_df.csv.gz',
    'monte_het_clean_clean_df.csv.gz',
    'monte_het_extend_100per_vpro_EXOENZYME_clean_df.csv.gz',
    'monte_het_extend_100per_vpro_MIXOTROPH_clean_df.csv.gz',
    'monte_het_extend_100per_vpro_OVERFLOW_clean_df.csv.gz',
    'monte_het_extend_100per_vpro_ROS_clean_df.csv.gz',
    'monte_het_extend_EXOENZYME_clean_df.csv.gz',
    'monte_het_extend_MIXOTROPH_clean_df.csv.gz',
    'monte_het_extend_OVERFLOW_clean_df.csv.gz',
    'monte_het_extend_ROS_clean_df.csv.gz',
    'monte_het_multi_clean_df.csv.gz',
    'monte_het_round2_100per_vpro_ROS_clean_df.csv.gz',
    'monte_ROS_round2_het_clean_df.csv.gz',
]


# In[ ]:


fnames_sum = [
    'monte_add_het_clean_clean_sum.csv.gz',
    'monte_het_add_100per_vpro_EXOENZYME_clean_sum.csv.gz',
    'monte_het_add_100per_vpro_OVERFLOW_clean_sum.csv.gz',
    'monte_het_add_100per_vpro_ROS_clean_sum.csv.gz',
    'monte_het_clean_clean_sum.csv.gz',
    'monte_het_extend_100per_vpro_EXOENZYME_clean_sum.csv.gz',
    'monte_het_extend_100per_vpro_MIXOTROPH_clean_sum.csv.gz',
    'monte_het_extend_100per_vpro_OVERFLOW_clean_sum.csv.gz',
    'monte_het_extend_100per_vpro_ROS_clean_sum.csv.gz',
    'monte_het_extend_EXOENZYME_clean_sum.csv.gz',
    'monte_het_extend_MIXOTROPH_clean_sum.csv.gz',
    'monte_het_extend_OVERFLOW_clean_sum.csv.gz',
    'monte_het_extend_ROS_clean_sum.csv.gz',
    'monte_het_multi_clean_sum.csv.gz',
    'monte_het_round2_100per_vpro_ROS_clean_sum.csv.gz',
    'monte_ROS_round2_het_clean_sum.csv.gz',
]
sum_df = pd.concat([pd.read_csv(os.path.join(dpath, f)) for f in fnames_sum], ignore_index=True)
sim_df = pd.concat([pd.read_csv(os.path.join(dpath, f)) for f in fnames_df], ignore_index=True)
sim_df = sim_df.drop(columns=['Unnamed: 0'])

minmse_df = pd.concat([pd.read_csv(os.path.join(dpath, f)) for f in fnames_predicted_classes], ignore_index=True)
minmse_df.VPRO = minmse_df.VPRO.str.replace('_monte_', '' ,regex=False)

# vpro_df = pd.read_csv('../ML_model/versatile_vpros.csv')
# vpro_df = vpro_df.loc[vpro_df.Versatile]

# minmse_df_vers = pd.merge(
#    vpro_df[['model', 'Phase', 'VPRO', 'Versatile',]],
#    minmse_df[['run_id', 'y_pred',  'VPRO', 'Sample', 'max_prob']],
#    on='VPRO', how='left')


df = pd.merge(minmse_df_vers, sim_df, on='run_id', how='left')


def gen_df_for_nas(y, sum_df=sum_df1, which_organism='all'):
    run_id = y.name
    param_vals = sum_df.loc[sum_df['run_id'].isin([run_id])].squeeze().to_dict()
    (var_names, init_var_vals, intermediate_names, calc_dydt, prepare_params_tuple
        ) = get_constants_per_organism(pro99_mode=False, which_organism=which_organism)
    par_tuple = prepare_params_tuple(param_vals)
    d = y.copy()
    d[intermediate_names] = d[var_names].apply(
        lambda x : calc_dydt(0, x.to_numpy(), par_tuple=par_tuple, return_intermediate=True), axis=1,
        result_type='expand')
    return d

def fill_NAs():
''' calculate values for params where we didn't print it (old versions) '''
    nacols = ['ROSproductionp', 'ROSproductionh', 'ROSlossp', 'ROSlossh',
           'deathbiomassNp', 'deathbiomassNh', 'deathstoreNp', 'deathstoreNh',
           'deathstoreCp', 'deathstoreCh', 'DON2DIN_exop', 'DON2DIN_exoh',
           'DON2DIN', 'additionalLossRatep', 'additionalLossRateh', 'regNp',
           'regNh', 'regCp', 'regCh', 'deathC_DOCp', 'deathC_DOCh', 'deathN_DONp',
           'deathN_DONh']


    df_of_nas = df.loc[df[nacols].isna().any(axis=1)]
    df_no_nas = df.loc[~df.run_id.isin(df_of_nas.run_id.unique())].copy()
    sum_df1 = sum_df.loc[sum_df.run_id.isin(df_of_nas.run_id.unique())].copy()


    df_filled = df_of_nas.groupby('run_id',group_keys=True).apply(gen_df_for_nas)
    df_final = pd.concat([df_filled, df_no_nas], ignore_index=True)
    return df_final
    


# # integrate

df = df.reset_index(drop=True)



from scipy import integrate
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
total_gross_uptake = df.groupby(['run_id', 'y_pred', 'model']).apply(lambda x: _integrate_column(x, collist))
total_gross_uptake = total_gross_uptake.reset_index()


N_uptake_cols = ['gross_uptakeINp', 'gross_uptakeINh', 'gross_uptakeONp', 'gross_uptakeONh']
total_gross_uptake['PP'] = total_gross_uptake['gross_uptakeICp']
total_gross_uptake['Total N uptake'] = total_gross_uptake[N_uptake_cols].sum(axis=1)


df['total C biomass'] = df['Bptotal[C]'] + df['Bhtotal[C]']
df['total fixed C'] = df['total C biomass'] + df['DOC'] + df['RDOC']


df = df.reset_index(drop=True)


last_day_df = df.loc[df.groupby('run_id').day.idxmax()]



id_keys = ['run_id', 'y_pred', 
           'model',
          ]
comb_df = pd.merge(
    last_day_df, 
    total_gross_uptake, 
    on=id_keys,
    suffixes=[' Final', ' Integrated'],
)


comb_df['PP / total fixed C'] = comb_df['PP'] / comb_df['total fixed C'] 
comb_df['N reuse'] = comb_df['Total N uptake'] / (INIT_DIN + INIT_DON)


comb_df.to_csv(os.path.join(dpath, 'integrated_last_core_models_versatile_het.csv.gz'), index=False)



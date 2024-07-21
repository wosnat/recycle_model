#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import math
import os

# import sys
# sys.path.append('..')
# from model_equations_separate_NC_store_numba import *


def _get_stable_range(x, col, startday=0):
    x1 = x.loc[x['day'].ge(startday)].copy()
    std_df = x1[col].rolling(window=3,center=True).std()
    #diff_df = tdf.QCp.diff()
    
    #Create a grouping series that tracks changes.
    groupr = std_df.ge(1e-2).cumsum()
    
    #Create a mapping dictionary that can translate from the groupr key to the actual value in the df.value column.
    #mapper = dict(zip(groupr, tdf.QCp))
    
    #Now we group and use ptp and nlargest. 
    #Finally, we use rename and mapper to translate the index value, which is the groupr value, back to the value value (phew, that's a tad confusing).
    #tdf.groupby(groupr).day.apply(np.ptp).nlargest(2).rename(mapper)
    range_df = x1.groupby(groupr).day.agg(['first', 'last', 'count'])
    range_df['range'] = range_df['last'] -  range_df['first'] 
    f, l, c, r = range_df.nlargest(1, 'range').squeeze()
    #print(range_df.nlargest(1, 'range'))
    #print (f,l,r)
    s,m = x1.loc[x1.day.between(f, l, inclusive='both'), col].agg(['std', 'mean'])
    
    return (f, l, m ,s, c, r)
    


def _get_stages(x):
    x1 = x.loc[x.day.gt(1)].reset_index(drop=True).copy()
    maxday = x1.day.max()
    maxPC = x1.loc[x1['Bptotal[C]'].idxmax(), 'day']
    maxHC = x1.loc[x1['Bhtotal[C]'].idxmax(), 'day'] 
    maxPN = x1.loc[x1['Bptotal[N]'].idxmax(), 'day'] 
    maxHN = x1.loc[x1['Bhtotal[N]'].idxmax(), 'day'] 
    minDIC = x1.loc[x1['DIC'].idxmin(), 'day'] 
    
    deadP = x1.loc[x1['Bptotal[C]'].le(2) & x1.day.gt(maxPC), 'day'].min()
    if deadP is np.NaN:
        deadP = maxday
    deadH = x1.loc[x1['Bhtotal[C]'].le(2)& x1.day.gt(maxHC), 'day'].min()
    if deadH is np.NaN:
        deadH = maxday
    first_QCp, last_QCp, mean_QCp ,std_QCp, count_QCp, length_QCp = _get_stable_range(x, 'QCp', 0)
    first_QCh, last_QCh, mean_QCh ,std_QCh, count_QCh, length_QCh  = _get_stable_range(x, 'QCh', 0)
    

    return pd.Series({
        'Max day Pro [C]' : maxPC, 
        'Max day Het [C]' : maxHC, 
        'Max day Pro [N]' : maxPN, 
        'Max day Het [N]' : maxHN, 
        'Min DIC day' : minDIC, 
        'Death day Pro' : deadP , 
        'Death day Het' : deadH ,
        
        'Equilibrium First day Pro' : first_QCp, 
        'Equilibrium last day Pro' : last_QCp, 
        'Equilibrium mean C/N Pro' : mean_QCp, 
        'Equilibrium std C/N Pro' : std_QCp, 
        'Equilibrium count Pro' : count_QCp, 
        'Equilibrium length Pro' : length_QCp, 
        
        'Equilibrium First day Het' : first_QCh, 
        'Equilibrium last day Het' : last_QCh, 
        'Equilibrium mean C/N Het' : mean_QCh, 
        'Equilibrium std C/N Het' : std_QCh, 
        'Equilibrium count Het' : count_QCh, 
        'Equilibrium length Het' : length_QCh, 
        
    })







def _mark_stages(x):
    #print (x.name)
    stages = stage_df.loc[stage_df.run_id.isin([x.name])].squeeze()
    #print(stages)
    #print([x.day.min()-1, stages['maxPC'], stages['stableQCP'], stages['deadP'], x.day.max()+1])
    maxPC = stages['Max day Pro [C]']
    pro_eq_count = stages['Equilibrium count Pro']
    pro_eq_start = stages['Equilibrium First day Pro']
    pro_eq_end = stages['Equilibrium last day Pro']
    het_eq_count = stages['Equilibrium count Het']
    eq_start = stages[['Equilibrium First day Pro', 'Equilibrium First day Het']].max()
    eq_end = stages[['Equilibrium last day Pro', 'Equilibrium last day Het']].min()

    # equilibrium only after maxPC
    eq_start = max(maxPC, eq_start)
    eq_end = max(maxPC, eq_end)
    pro_eq_start = max(maxPC, pro_eq_start)
    pro_eq_end = max(maxPC, pro_eq_end)
    
    # the defualt set of stages, assuming 2 equilibrium states. (what abount equilibrium late pro?)    
    bins = [x.min()-1, maxPC, pro_eq_start, eq_start, eq_end, pro_eq_end, x.max()+1]
    labels=['Bloom', 'Bust', 'Equilibrium Pro', 'Equilibrium', 'Equilibrium Pro (cont)', 'Post equilibrium']
    if (pro_eq_count <= 5):
        # no pro equilibrium
        bins = [x.min()-1, maxPC,  x.max()+1]
        labels=['Bloom', 'Bust', ]
    elif (eq_start + 5>= eq_end) or (het_eq_count <= 5):
        # no common equilibrium
        bins = [x.min()-1, maxPC, pro_eq_start, pro_eq_end, x.max()+1]
        labels=['Bloom', 'Bust',  'Equilibrium Pro', 'Post equilibrium']
    elif pro_eq_start >= eq_start:
        # no early Pro equilibrium 
        bins = [x.min()-1, maxPC, eq_start, eq_end, pro_eq_end, x.max()+1]
        labels=['Bloom', 'Bust',  'Equilibrium', 'Equilibrium Pro (cont)', 'Post equilibrium']
    elif pro_eq_end <= eq_end:
        bins = [x.min()-1, maxPC, eq_start, eq_end, x.max()+1]
        labels=['Bloom', 'Bust',  'Equilibrium', 'Post equilibrium']
        
        
    labels = [labels[i] for i in range(len(labels)) if bins[i] != bins[i+1]]
    #print(bins)
    #print(labels)
    return  pd.cut(
        x,
        bins=bins,
        labels=labels,
        ordered=False, duplicates='drop',
    )



if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='cleanup PONLY runs and produce VPROS.')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--infpath", help="input dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    
    
    args = parser.parse_args()
    out_dpath = args.outdpath
    if out_dpath != '':
        os.makedirs(out_dpath, exist_ok=True)

#######################

    fpath = args.infpath
    # results/final/het/monte_het_add_100per_vpro_ROS_clean_df.csv.gz
    fname = os.path.basename(fpath)
    session_id = fname.replace('_df.csv.gz', '')
    
    session_id = args.run_id
    df = pd.read_csv(fpath)
    df = df.reset_index(drop=True)


    per_sec_cols = [
         'gross_uptakeINp',
           'gross_uptakeINh', 'gross_uptakeONp', 'gross_uptakeONh',
           'gross_uptakeICp', 'gross_uptakeICh', 'gross_uptakeOCp',
           'gross_uptakeOCh', 'uptakeNp', 'uptakeNh', 'uptakeCp', 'uptakeCh',
            'biosynthesisNp', 'biosynthesisNh', 'respirationCp',
           'respirationCh', 'biomass_breakdownCp', 'biomass_breakdownCh',
           'overflowNp', 'overflowNh', 'overflowCp', 'overflowCh', 
           
           'ROSproductionp', 'ROSproductionh', 'ROSlossp', 'ROSlossh',
           'deathbiomassNp', 'deathbiomassNh', 'deathstoreNp', 'deathstoreNh',
           'deathstoreCp', 'deathstoreCh', 'DON2DIN_exop', 'DON2DIN_exoh',
           'DON2DIN', 'additionalLossRatep', 'additionalLossRateh',  'deathC_DOCp', 'deathC_DOCh', 'deathN_DONp',
           'deathN_DONh'
    ]
    for c in per_sec_cols:
        df[c] = df[c]*seconds_in_day



    df['alive'] = 'alive'
    maskp = df['Bptotal[C]'].le(2)
    maskh = df['Bhtotal[C]'].le(2)
    df.loc[maskp, 'alive'] =  'deadp'
    df.loc[maskh, 'alive'] =  'deadh'
    df.loc[maskp & maskh, 'alive'] =  'dead'


    df['Het/Pro C'] = df['Bhtotal[C]'].div(df['Bptotal[C]'])
    df['Het/Pro N'] = df['Bhtotal[N]'].div(df['Bptotal[N]'])

    stage_df = df.groupby(['run_id',]).apply(_get_stages).reset_index()

    df['stage'] = df.groupby('run_id')['day'].transform(_mark_stages)


    stat_df = df.groupby(['run_id', 'stage']).agg({
        'DON' : ['mean', 'std'],  
        'DIN' : ['mean', 'std'],
           'DOC': ['mean', 'std'],
        'ROS': ['mean', 'std'],
        'QCp': ['mean', 'std'],
        'QCh': ['mean', 'std'],
        'Het/Pro C': ['mean', 'std'],
        'Het/Pro N': ['mean', 'std'],
        'day': lambda x : x.max() - x.min(),
        'alive' : lambda x : ' '.join(x.unique()),
    })

    stat_df.columns = stat_df.columns.get_level_values(0) + ' ' +  stat_df.columns.get_level_values(1) 
    stat_df.rename(columns={'day <lambda>' : 'length (days)', 'alive <lambda>' : 'alive states'}, inplace=True)


    pstat_df = stat_df.reset_index().pivot(
        index=['run_id'],
        columns='stage',
    )
    pstat_df.columns = pstat_df.columns.get_level_values(0) + ' [' +  pstat_df.columns.get_level_values(1) + ']'
    pstat_df = pstat_df.reset_index()


    merged_stage_df = pd.merge(stage_df, pstat_df, on=['run_id'],)
    merged_stage_df

    merged_stage_df.to_csv(os.path.join(out_dpath, f'stage_stats_{session_id}.csv'))

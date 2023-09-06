#!/usr/bin/env python
# coding: utf-8

import os
import glob
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
from scipy.integrate import odeint
from sklearn.metrics import mean_squared_error
import json
import subprocess
import sys
import re
import time

import calc_csat

from numba import jit



# Redfield ratio
R_CN = 6.625

# parameter values
pro_radius = 0.3628;  # "MED4" = 9312
pro_radius = 0.5;  # "MED4" = 9312
alt_radius = 1;  # "MED4" = 9312
pro_vol = (4 / 3) * math.pi * pro_radius ** 3;
alt_vol = (4 / 3) * math.pi * alt_radius ** 3;
pro_alt_vol_ratio = alt_vol / pro_vol
seconds_in_day = 24*60*60

# umol/cell
Qmaxp = 1.5e-9
Qminp = 7e-10
Qmaxh = 1.5e-9 * pro_alt_vol_ratio
Qminh = 7e-10  * pro_alt_vol_ratio
Qh = (Qminh + Qmaxh) / 2
Qp = (Qminp + Qmaxp) / 2

# H params
# H: C: 30 - 292 fg cell-1
# H: N: 25-50 mean 37 fg cell-1
# H C:N 4-6.6 mean 4.9 fg cell-1


# P params
# N 5-7 fg/cell
# C 20-50 fg/cell
# C:N 5.5-7, 5-12
# michal 10-21 N
# shira/yara 17-30 N
# Dalit N 11.8-22.3 fg cell-1 N
# Dalit C 71-134 fg cell-1
# Dalit C:N 6

# parameter values
R_P = 7
R_H = 4.5 # in the range 4-6 (for DSS N limited, 11)

# fg -> umol 14 (N mulecular weight) * 1e-9 (fmol -> umol)
Qp = 25  * 1e-9 / 14
Qh = 40 * 1e-9 / 14



# DIC exchange
h = 0.115 # height in m of media in the tube
Kg = 1.40E-06; #% m sec-1 measured in the lab
# Kg = 0.4968 / seconds_in_day # m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
#c_sat = INIT_DIC

# from Michal:
# % CO2 parameters
# Kg = 0.4968/86400; #% m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998

B  = 10;           #% Revelle buffer factor, Mick's book
 
rhoref = 1024.5;        #% reference density of seawater
thetaK = 273.15 + 22.0; #% temperature (degrees K) 
salt   = 38.5;          #% salinity (PSU) 
pt     = 50 * rhoref / 1e6;  #% inorganic phosphate (mol/^3) 
sit    = 0.0;           #% inorganic silicate (mol/^3) 
pCO2   = 400e-6;        # atmospheric reference pCO2 level (atmospheres)
                        # for which to find equilibrium dic, csat
ta      = 2650 * rhoref / 1e6; # total alkalinity (eq/m^3), from Dalit
# (alk from To new PRO99 medium. addition of 1mM BC)

# equilibrium total inorganic carbon (mol/m^3)
csatd = calc_csat.calc_csat(T=thetaK, S=salt, pco2eq=pCO2, pt=pt, sit=sit, ta=ta)
c_sat = csatd * 1e6 / rhoref; 
#    h = (.09/-27/86400)*i+.105; % the hight of the water columns in meters
#    
#    % alkalinity and Csat calculation
#    ta1     = (112.46/86400)*i + 2617.6;
#    ta      = ta1 * rhoref / 1e6;
#    [csatd] = calc_csat(thetaK, salt, pCO2, pt, sit, ta);
#    Csat    = csatd * 1e6 / rhoref;

# initial concentrations
INIT_DIN = 100
INIT_DIN_PRO99 = 800
INIT_DON = 20
INIT_RDON = 0
INIT_RDOC = 0
# Dalit DIC: 1618.825333  or 1.62E+03 uM
INIT_DIC = c_sat

# Dalit: DOC 0.047155376888899 nM ?
# Dalit init TOC 16 mM
INIT_DOC = INIT_DON * R_CN
INIT_BP = 1e9 * Qp
INIT_BH = 1e10 * Qh
INIT_CP = 0
INIT_CH = 0
INIT_NP = 0
INIT_NH = 0

INIT_BH_CC = 5e9 * Qh # the actual concentration in the measurements
INIT_ROS = 0.2 # Morris, J. Jeffrey, et al. "Dependence of the cyanobacterium Prochlorococcus on hydrogen peroxide scavenging microbes for growth at the ocean's surface." PloS one 6.2 (2011): e16805.
INIT_SP = 0
INIT_SH = 0

def get_param_vals(model_name):
    param_df = pd.read_excel( 'Model_Parameters Store model.xlsx',)
    param_df['values'] = param_df['full model']
    if model_name not in [None, 'FULL']:
        model_col = f'{model_name} model hardcoded parameters'
        param_df.loc[~param_df[model_col].isna(), 'values'] = param_df.loc[~param_df[model_col].isna(), model_col] 
    param_vals_series = param_df['values']
    param_vals_series.index = param_df.parameter
    return param_vals_series.to_dict()
    

    # K’s (affinity).
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake


param_vals = get_param_vals(None)

def get_param_tuning_values(model_name, organism_to_tune):
    """ return 
    params_to_update (list of strings)
    bounds (list of pairs)
    log_params (list of bools
    """
    param_df = pd.read_excel( 'Model_Parameters.xlsx',)
    org_col = f'Tunable parameters ({organism_to_tune} fitting)'

    if model_name not in [None, 'FULL']:
        # restrict to given model
        model_col = f'{model_name} model hardcoded parameters'
        tunable_param_df = param_df.loc[param_df[model_col].isna() & param_df[org_col].isin(['Yes'])] 
    else: 
        # take the full model
        tunable_param_df = param_df.loc[param_df[org_col].isin(['Yes'])] 

    params_to_update = list(tunable_param_df['parameter'])    
    bounds = list(zip(tunable_param_df['lower bound'],tunable_param_df['upper bound']))    
    
    log_params = list(tunable_param_df['logscale fitting'].map({'Yes': True, 'No': False})) 

    return (params_to_update, bounds, log_params)

    

    # K’s (affinity).
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake


param_vals = get_param_vals(None)
     

def print_params(param_vals=param_vals):
    for i in param_vals:
        print(i, f' = {param_vals[i]:.2e}, {param_vals[i] * seconds_in_day:.2e}')

def basic_model_ode(time, var_vals, param_vals):
    # parameters
    pars = param_vals
    paramCN  = np.array([pars['Rp'], pars['Rh']])
    paramQCmax  = np.array([pars['QCmaxp'], pars['QCmaxh']])
    paramQCmin  = np.array([pars['QCminp'], pars['QCminh']])
    paramKmtb  = np.array([pars['Kmtbp'], pars['Kmtbh']])
    paramb  = np.array([pars['bp'], pars['bh']])
    paramr0  = np.array([pars['r0p'], pars['r0h']])
    paramM = np.array([pars['Mp'], pars['Mh']])
    paramE_leak = np.array([pars['E_leakp'], pars['E_leakh']])
    paramE_ROS = np.array([pars['E_ROSp'], pars['E_ROSp']])
    paramVmaxROSh = pars['VmaxROSh']
    paramK_ROSh = pars['K_ROSh']
    paramgamma_DON2DIN = np.array([pars['gamma_DON2DINp'], pars['gamma_DON2DINh']])
    paramgammaD = np.array([pars['gammaDp'], pars['gammaDh']])
    paramROS_decay = pars['ROS_decay']
    kns = np.array(
        [[pars['KINp'], pars['KONp'], pars['KICp'], pars['KOCp']],
        [pars['KINh'], pars['KONh'], pars['KICh'], pars['KOCh'], ]]
    ).T
    vmax = np.array(
        [[pars['VmaxINp'], pars['VmaxONp'], pars['VmaxICp'], pars['VmaxOCp']],
        [pars['VmaxINh'], pars['VmaxONh'], pars['VmaxICh'], pars['VmaxOCh'], ]]
    ).T

    # variables
    var_names  = ['Bp', 'Np',  'Cp',  'Bh',  'Nh',  'Ch',   'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    vars = {n:v for n,v in zip(var_names,var_vals)}
    DIC = vars['DIC']
    Bp = vars['Bp']
    Np = vars['Np']
    Cp = vars['Cp']
    Bh = vars['Bh']
    Nh = vars['Nh']
    Ch = vars['Ch']
    DIN = vars['DIN']
    DON = vars['DON']
    RDON = vars['RDON']
    DIC = vars['DIC']
    DOC = vars['DOC']
    RDOC = vars['RDOC']
    ROS = vars['ROS']
    
    resources = np.array(
        [DIN, DON, DIC, DOC, ]
    ).reshape(4, 1)
    DIN_IDX = 0
    DON_IDX = 1
    DIC_IDX = 2
    DOC_IDX = 3
    

    biomass = np.array([Bp, Bh])
    storeN  = np.array([Np, Nh])
    storeC  = np.array([Cp, Ch])

    ROSdecay = ROS * paramROS_decay
    netROS = ROS - ROSdecay

    QC = (storeC + biomass * paramCN) / (storeN + biomass)
  
    # monod ratios
    limits = resources / (resources + kns)
    # regulate uptake by not letting the stores grow too large 
    # question: is DIC uptake by photosynthesis regulated?
    regC = 1 - (QC / paramQCmax)
    regN = 1 - (paramQCmin / QC)
    reg =  np.vstack([
        np.atleast_2d(regN).repeat(repeats=2, axis=0), 
        np.atleast_2d(regC).repeat(repeats=2, axis=0)
    ])
    reg_clipped = np.clip(reg, a_min =0.0, a_max=1.0)
    # gross uptake (regardless of C:N ratio)
    # vmax = muinfp* VmaxIp / Qp
    # umol N /L  or umol C /L  
    gross_growth = biomass * vmax * limits  * reg_clipped 
    uptakeN = gross_growth[DIN_IDX,:] + gross_growth[DON_IDX,: ],   
    uptakeC = gross_growth[DIC_IDX,:] + gross_growth[DOC_IDX,: ],   

    # net uptake (maintains C:N ratios)
    # umol N / L
    biosynthesisN = np.minimum(storeN + uptakeN, (storeC + uptakeC) / paramCN)* paramKmtb
    # Respiration – growth associated bp/bh and maintenance associated r0p/r0h
    # b * growth + r0 * biomass
    # umol C/L
    respirationC = (paramb * biosynthesisN + biomass * paramr0) * paramCN
    # if C store is not big enough, break some of the biomass into the stores
    # make sure Cp is not negative
    # umol C/L
    biomass_breakdownC = np.clip(
        respirationC +  biosynthesisN * paramCN - storeC - uptakeC, 
        a_min=0, a_max=None,
    )

    # store change - uptake minus biosynthesis and respiration
    netDeltaN = uptakeN + biomass_breakdownC / paramCN - biosynthesisN 
    netDeltaC = uptakeC + biomass_breakdownC - biosynthesisN * paramCN  - respirationC

    # overflow -
    # make the store maintain the C:N ratio and exude the rest
    # Overflow quantity 
    # Oh/Op: enable overflow (0 or 1)
    # umol N / L
    store_keepN = np.minimum(storeN + netDeltaN, (storeC + netDeltaC) / paramCN) 
    overflowN = storeN + netDeltaN - store_keepN
    overflowC = storeC + netDeltaC - store_keepN * paramCN

    dic_air_water_exchange   = - (DIC - c_sat) / ((h * DIC) / (Kg * B * 0.01 * DIC))
    # death
    # Need to explain why we used exponential decay – in ISMEJ we show that other formulations are better for co-cultures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.
    deathN = paramM * biomass
    # leakiness formulated as fraction of biomass (“property tax”)
    # We assume that the vast majority of C and N biomass is in organic form, hence leakiness is to organic. We assume that overflow is also to organic in both organisms, as for the phototroph this is the release of fixed C (or inorganic N incorporated into e.g. AA) which cannot be used for growth. For the heterotrophs we assume overflow metabolism to be the inefficient use of organic C (e.g. not fully oxidized) to maximize growth rate (*citation E coli).
    leakinessN = paramE_leak * biomass
    # ROS production depends on biomass
    ROSrelease = paramE_ROS * biomass
    ROSbreakdownh = VmaxROSh * netROS / (netROS + K_ROSh) * Bh
    # DIN breakdown due to exoenzymes
    DON2DIN = paramgamma_DON2DIN * biomass * DON
    # final differential equations
    dBdt = biosynthesisN - biomass_breakdownC / paramCN - deathN - leakinessN 
    dNdt = netDeltaN - overflowN
    dCdt = netDeltaC - overflowC
    dDONdt = np.sum(
        deathN * paramgammaD + leakinessN - gross_growth[DON_IDX,:] - DON2DIN
    )
    dDOCdt = np.sum(
        deathN * paramgammaD * paramCN + 
        leakinessN * paramCN + overflowC - gross_growth[DOC_IDX,:])

    # In discussion can state that if DIN is produced also through overflow or leakiness then this could support Pro growth, but this is not encoded into our model.
    # Assuming that recalcitrant DON isd released only during mortality (discuss release through leakiness)
    # Assuming RDON/RDOC is recalcitrant to both organisms    dRDONdt = sum(deathN * (1 - paramgammaD))
    dRDOCdt = sum(deathN * paramCN * (1 - paramgammaD))
    dDINdt = sum(overflowN + DON2DIN - gross_growth[DIN_IDX,:])
    dDICdt = dic_air_water_exchange + sum(respirationC - gross_growth[DIC_IDX,:])
    dROSdt = sum(ROSrelease) - ROSdecay - ROSbreakdownh
    return ([dBdt[0], dNdt[0], dCdt[0], dBdt[1], dNdt[1], dCdt[1], dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt])


# functions

# number of cells (Q is a parameter)
# cells/L
#Xp = Bp / Qp
#Xh = Bh / Qh


def print_equations():
    sfunc_list = [dBpdt, dNpdt, dCpdt, dBhdt, dNhdt, dChdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dABpdt, dABhdt]
    var_names  = ['Bp', 'Np',  'Cp',  'Bh',  'Nh',  'Ch',   'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS',   'ABp',  'ABh']
    for n,f in zip(var_names, sfunc_list):
        print(f'd{n}/dt')
        print(str(f))

def get_main_init_vars(pro99_mode):
    if pro99_mode:        
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_BH,INIT_NH,INIT_CH,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    else:
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_BH,INIT_NH,INIT_CH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    return init_vars

def get_main_data(param_vals_str, pro99_mode):
    sfunc_list = [dBpdt, dNpdt, dCpdt, dBhdt, dNhdt, dChdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dABpdt, dABhdt]
    var_list   = [ Bp,   Np,    Cp,    Bh,    Nh,    Ch,     DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS,    ABp,    ABh]
    var_names  = ['Bp', 'Np',  'Cp',  'Bh',  'Nh',  'Ch',   'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS',   'ABp',  'ABh']
    init_vars = get_main_init_vars(pro99_mode)
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    #final_func = lambdify(var_list, subs_funclist, modules=['math'])
    #final_func_jit = jit(final_func, nopython=True)  
    final_func_jit = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func_jit(*y)

    interm_sfunc_list = [
        gross_uptakeINp, 
        gross_uptakeONp, 
        gross_uptakeINh, 
        gross_uptakeONh, 
        gross_uptakeICp, 
        gross_uptakeOCp, 
        gross_uptakeICh, 
        gross_uptakeOCh, 
        uptakeNp, 
        uptakeNh, 
        uptakeCp, 
        uptakeCh, 
        regQCp,
        regQCh,
        regQNp,
        regQNh,
        bio_synthesisN_p, 
        bio_synthesisN_h, 
        respirationCp, 
        respirationCh, 
        biomass_breakdown_for_respirationCp, 
        biomass_breakdown_for_respirationCh, 
        netDeltaNp, 
        netDeltaNh, 
        netDeltaCp, 
        netDeltaCh, 
        store_keepNp,
        store_keepNh,
        overflowNp, 
        overflowNh, 
        overflowCp, 
        overflowCh, 
        dic_air_water_exchange,
        ABreleasep,
        ABreleaseh, 
        death_ratep, 
        death_rateh, 
        deathp, 
        deathh, 
        leakinessOp, 
        leakinessIp, 
        leakinessOh, 
        leakinessIh, 
        ROSreleasep,
        ROSreleaseh, 
        ROSbreakdownh, 
        DON2DINp, 
        DON2DINh, 
    ]
    interm_names = [
        'gross_uptakeINp', 
        'gross_uptakeONp', 
        'gross_uptakeINh', 
        'gross_uptakeONh', 
        'gross_uptakeICp', 
        'gross_uptakeOCp', 
        'gross_uptakeICh', 
        'gross_uptakeOCh', 
        'uptakeNp', 
        'uptakeNh', 
        'uptakeCp', 
        'uptakeCh', 
        'regQCp',
        'regQCh',
        'regQNp',
        'regQNh',
        'bio_synthesisN_p', 
        'bio_synthesisN_h', 
        'respirationCp', 
        'respirationCh', 
        'biomass_breakdown_for_respirationCp', 
        'biomass_breakdown_for_respirationCh', 
        'netDeltaNp', 
        'netDeltaNh', 
        'netDeltaCp', 
        'netDeltaCh', 
        'store_keepNp',
        'store_keepNh',
        'overflowNp', 
        'overflowNh', 
        'overflowCp', 
        'overflowCh', 
        'dic_air_water_exchange',
        'ABreleasep',
        'ABreleaseh', 
        'death_ratep', 
        'death_rateh', 
        'deathp', 
        'deathh', 
        'leakinessOp', 
        'leakinessIp', 
        'leakinessOh', 
        'leakinessIh', 
        'ROSreleasep',
        'ROSreleaseh', 
        'ROSbreakdownh', 
        'DON2DINp', 
        'DON2DINh', 
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    #intermediate_func = lambdify(var_list, interm_funclist, modules=['math'])
    #intermediate_func = jit(intermediate_func, nopython=True)  
    intermediate_func = lambdify(var_list, interm_funclist)
    #interm_names = []
    #intermediate_func = None
    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_init_vars(pro99_mode):
    if pro99_mode:        
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    else:
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    return init_vars
    
def get_ponly_data(param_vals_str, pro99_mode):
    sfunc_list = [dBpdt,  dNpdt, dCpdt, dDONdt_ponly, dRDONdt_ponly, dDINdt_ponly, dDOCdt_ponly, dRDOCdt_ponly, dDICdt_ponly, dROSdt_ponly, dABpdt, dABhdt_ponly]
    var_list   = [ Bp,    Np,    Cp,    DON,    RDON,    DIN,    DOC,    RDOC,  DIC,    ROS,    ABp,   ABh]
    var_names  = ['Bp',  'Np',  'Cp',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'ABp', 'ABh']
    init_vars = get_ponly_init_vars(pro99_mode)
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    #final_func = lambdify(var_list, subs_funclist, modules=['math'])
    #final_func_jit = jit(final_func, nopython=True)  
    final_func_jit = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func_jit(*y)

    interm_sfunc_list = [
        gross_uptakeINp, 
        gross_uptakeONp, 
        gross_uptakeICp, 
        gross_uptakeOCp, 
        uptakeNp, 
        uptakeCp, 
        regQCp,
        regQNp,
        bio_synthesisN_p, 
        respirationCp, 
        biomass_breakdown_for_respirationCp, 
        netDeltaNp, 
        netDeltaCp, 
        store_keepNp,
        overflowNp, 
        overflowCp, 
        dic_air_water_exchange,
        ABreleasep,
        death_ratep, 
        deathp, 
        leakinessOp, 
        leakinessIp, 
        ROSreleasep,
        DON2DINp, 
    ]
    interm_names = [
        'gross_uptakeINp', 
        'gross_uptakeONp', 
        'gross_uptakeICp', 
        'gross_uptakeOCp', 
        'uptakeNp', 
        'uptakeCp', 
        'regQCp',
        'regQNp',
        'bio_synthesisN_p', 
        'respirationCp', 
        'biomass_breakdown_for_respirationCp', 
        'netDeltaNp', 
        'netDeltaCp', 
        'store_keepNp',
        'overflowNp', 
        'overflowCp', 
        'dic_air_water_exchange',
        'ABreleasep',
        'death_ratep', 
        'deathp', 
        'leakinessOp', 
        'leakinessIp', 
        'ROSreleasep',
        'DON2DINp', 
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    #intermediate_func = lambdify(var_list, interm_funclist, modules=['math'])
    #intermediate_func = jit(intermediate_func, nopython=True)  
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_init_vars(pro99_mode):
    if pro99_mode:        
        init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    else:
        init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    return init_vars
    
def get_honly_data(param_vals_str, pro99_mode):

    sfunc_list = [dBhdt, dNhdt, dChdt, dDONdt_honly, dRDONdt_honly, dDINdt_honly, dDOCdt_honly, dRDOCdt_honly,dDICdt_honly, dROSdt_honly, dABpdt_honly, dABhdt]
    var_list   = [Bh,    DON,    Nh,    Ch,     RDON,    DIN,    DOC,    RDOC,   DIC,    ROS,    ABp,   ABh]
    var_names  = ['Bh',  'DON', 'Nh',  'Ch',   'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'ABp', 'ABh']
    init_vars = get_honly_init_vars(pro99_mode)
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist, modules=['math'])
    final_func_jit = jit(final_func, nopython=True)  
    calc_dydt = lambda t, y : final_func_jit(*y)

    interm_sfunc_list = [
    ]
    interm_names = [
    ]

    #interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    #intermediate_func = lambdify(var_list, interm_funclist, modules=['math'])
    #intermediate_func = jit(intermediate_func, nopython=True)  
    intermediate_func = None
    return var_names, init_vars, calc_dydt, interm_names, intermediate_func


def print_dydt0(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    for i,j, k in zip(var_names, dydt0, init_vars):
        print(f'd{i}/dt = {j:.2e}, init {i} = {k:.2e}, newval = {k+j:.2e}')


def print_intermediate0(intermediate_func, interm_names, init_vars):
    for i,j in zip(interm_names, intermediate_func(*init_vars)):
        print(f'{i:<4} = {j:.2e}')


def biomass_diff0(calc_dydt, var_names, init_vars, param_vals):
    dydt0 = calc_dydt(0, init_vars)
    R_P = param_vals['Rp']
    V = dict(zip(var_names, dydt0))
    print (f"dBp/dt + dBh/dt + dDON/dt + dRDON/dt + dDIN/dt = { V['Bp'] + V['Np'] + V['Bh'] + V['Nh'] + V['DON'] + V['RDON'] + V['DIN'] }")
    print (f"dBp/dt + dBh/dt + dDOC/dt + dRDOC/dt + dDIC/dt = { V['Bp']*R_P + V['Cp'] + V['Bh']*R_H + V['Ch'] + V['DOC'] + V['RDOC'] + V['DIC'] }")

def biomass_diff0_ponly(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    V = dict(zip(var_names, dydt0))
    print (f"dBp/dt  + dDON/dt + dRDON/dt + dDIN/dt = { V['Bp']  + V['DON'] + V['RDON'] + V['DIN'] }")

def biomass_diff0_honly(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    V = dict(zip(var_names, dydt0))
    print (f"dBh/dt + dDON/dt + dRDON/dt + dDIN/dt = {  V['Bh'] + V['DON'] + V['RDON'] + V['DIN'] }")



def run_solver_ivp(calc_dydt, init_vars, days, t_eval):
    tstart = 0
    if t_eval is None: 
        tend = days*seconds_in_day
        t_eval = np.arange(tstart, tend, 3600*4)
    else:
        tend= np.max(t_eval)
    sol = solve_ivp(
        fun=calc_dydt, y0=init_vars,
        t_span=[tstart, tend], t_eval=t_eval, max_step=100, first_step=1)
    #print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval={t_eval})')
    #print(sol.message)
    return sol


def solver2df_ivp(sol, var_names, interm_names, intermediate_func, param_vals):
    d = dict(zip(var_names, sol.y))
    d['t'] = sol.t
    df = pd.DataFrame(data=d)
    df['day'] = df['t']/seconds_in_day
    df[interm_names] = df[var_names].apply(lambda x : intermediate_func(*x), axis=1, 
                                           result_type='expand')
    if 'Bp' in df.columns:
        df['Bp[C]'] = df['Bp']*param_vals[str(Rp)]
    if 'Bh' in df.columns:
        df['Bh[C]'] = df['Bh']*param_vals[str(Rh)]
    if 'ABp' in df.columns:
        df['ABp[C]'] = df['ABp']*param_vals[str(Rp)]
    if 'ABh' in df.columns:
        df['ABh[C]'] = df['ABh']*param_vals[str(Rh)]
    return df


def run_solver_ode(calc_dydt, init_vars, days, t_eval):
    if t_eval is None: 
        tstart = 0
        tend = days*seconds_in_day
        t_eval = np.arange(tstart, tend, 3600*4)

        #print(t_eval)
    y = odeint(
        func=calc_dydt, y0=init_vars,
        t=t_eval, hmax=100, h0=1, tfirst=True)
    #print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval={t_eval})')
    #print(sol.message)
    return t_eval, y
    #return solver2df_ode(sol, t_eval, var_names, interm_names, intermediate_func)

def solver2df_ode(sol, var_names, interm_names, intermediate_func, param_vals, t_eval):
    d = dict(zip(var_names, sol[1].T))
    d['t'] = sol[0]
    df = pd.DataFrame(data=d)
    if t_eval is not None:
        df = df.loc[np.rint(df.t).isin(np.rint(t_eval))]
    df['day'] = df['t']/seconds_in_day
    if (interm_names):
        df[interm_names] = df[var_names].apply(lambda x : intermediate_func(*x), axis=1, 
                                           result_type='expand')
    if 'Bp' in df.columns:
        df['Bp[C]'] = df['Bp']*param_vals[str(Rp)]
    if 'Bh' in df.columns:
        df['Bh[C]'] = df['Bh']*param_vals[str(Rh)]
    if 'ABp' in df.columns:
        df['ABp[C]'] = df['ABp']*param_vals[str(Rp)]
    if 'ABh' in df.columns:
        df['ABh[C]'] = df['ABh']*param_vals[str(Rh)]
    return df


def run_solver(calc_dydt, init_vars, days=140, t_eval=None):
    tstart = time.process_time()
    sol = run_solver_ode(calc_dydt, init_vars, days, t_eval=t_eval)
    tend =  time.process_time()
    print ('simulation time', tend - tstart)
    return sol
    
def solver2df(sol, var_names, interm_names, intermediate_func, param_vals, t_eval=None):
    return solver2df_ode(sol, var_names, interm_names, intermediate_func, param_vals, t_eval)

def get_t_eval(maxday, step = 1000, ref_times = None):
    tstart = 0
    tend = maxday*seconds_in_day
    sim_times = np.arange(tstart, tend, step)
    if ref_times is not None:
        ref_times = np.rint(ref_times)
        sim_times = np.concatenate((sim_times, ref_times))
        sim_times = np.unique(sim_times)
    return sim_times

##################################################################
##################################################################
##################################################################

# run the model as a separate process to allow timeout
def params2json(params, fpath):
    with open(fpath, mode='w') as fp:
        json.dump(params,
                  fp, sort_keys=True, indent=4)

def json2params(default_params, fpath):
    with open(fpath) as fp:
        tparams = json.load(fp)
    new_param_vals = default_params.copy()
    new_param_vals.update(tparams)
    return new_param_vals


#def _rmse(df, refdf, refcol, col):
#    smallrefdf = refdf.dropna(subset=[refcol])
#    ref_t = np.rint(smallrefdf['t'])
#    tdf = df.loc[df.t.isin(ref_t)]
#    return mean_squared_error(tdf[col], smallrefdf[refcol])
   
def _mse(x, refdf, refcol, col, timecol, tolerance):
    #print(x.columns)
    tdf = pd.merge_asof(x[[timecol, col]].dropna(), refdf[[timecol, refcol]], on=timecol, direction='nearest', tolerance=tolerance).dropna()
    return pd.Series({
        'compare_points': tdf.shape[0], 
        'MSE': mean_squared_error(tdf[col], tdf[refcol])}
    )
def compute_mse(df, refdf, refcol= 'ref_Bp', col='Bp', timecol='t', tolerance=100):
    #print (x)
    mse_df = refdf.groupby(['Sample', 'full name', 'Group',]
                        ).apply(lambda y : _mse(df,y, refcol= refcol, col=col, timecol=timecol, tolerance=tolerance))
    mse_df = mse_df.reset_index()
    return mse_df

# this is sensitivity specific
# def compute_mse(df, refdf, refcol= 'ref_Bp', col='Bp', timecol='t', tolerance=100):
# return df.groupby(['sen_param', 'model', 'idx', 'run_id' ]
# ).apply(lambda x : _mse_all(x, refdf, refcol= refcol, col=col, timecol=timecol, tolerance=tolerance))    


def run_with_params_json(json_fpath_list, days, refdf, out_dpath, out_fprefix, which_organism, pro99_mode, t_eval):
    perr = -1
    orig_t_eval = t_eval
    new_params = param_vals
    for json_fpath in json_fpath_list:
        new_params = json2params(new_params, json_fpath)
    if which_organism == 'ponly':
        var_names, init_vars, calc_dydt, interm_names, intermediate_func = get_ponly_data(param_vals_str=new_params, pro99_mode=pro99_mode)
    elif which_organism == 'honly':
        var_names, init_vars, calc_dydt, interm_names, intermediate_func = get_honly_data(param_vals_str=new_params, pro99_mode=pro99_mode)
    else:
        var_names, init_vars, calc_dydt, interm_names, intermediate_func = get_main_data(param_vals_str=new_params, pro99_mode=pro99_mode)
    if t_eval is None:
        if refdf is not None:
            orig_t_eval = np.rint(refdf['t'].drop_duplicates().sort_values()).values
            t_eval = get_t_eval(days, ref_times = orig_t_eval)
        else:
            t_eval = get_t_eval(days)
    else:
        t_eval = get_t_eval(days, ref_times = t_eval)

    
    sol = run_solver(calc_dydt, init_vars, t_eval=t_eval, days=days)
    sumdf = pd.DataFrame({str(k): v for k,v in new_params.items()}, index=[0])
    sumdf['run_id'] = out_fprefix
    #sumdf['status'] = sol.status
    #if sol.status != 0:
    #    sumdf['message'] = sol.message
    #if sol.success:
    df = solver2df(sol, var_names, None, intermediate_func, new_params, t_eval=orig_t_eval)
    df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'), compression='gzip')

    if refdf is not None:
        try:
            mse_df = compute_mse(df, refdf, refcol= 'ref_Bp', col='Bp', timecol='t', tolerance=100)
            mse_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_mse.csv.gz'), compression='gzip')
            perr = mse_df['MSE'].min()

        except Exception as inst:
            print(inst)
            pass
    sumdf.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'), compression='gzip')
    return perr 
   
def generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60, which_organism='all', pro99_mode=False, t_eval=None):
    if pro99_mode:
        out_fprefix = f'{out_fprefix}_pro99'
    
    hash_val = str(hash(tuple(params.values())))
    run_id = f'{out_fprefix}_h{hash_val}'
    
    json_fpath = os.path.join(json_dpath, f'{run_id}_params.json')
    params2json(params, json_fpath)
    return run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout, which_organism, pro99_mode, t_eval)


def get_params(X, params_to_update, param_vals, log_params=None): 
    new_param_vals = param_vals.copy()
    if log_params is not None:
        X = [np.exp(x) if lg else x for x,lg in zip(X, log_params)]
    
    new_param_vals.update({k : v for k,v in zip(params_to_update, X)})
    return new_param_vals

def generate_json_and_run_from_X(X, params_to_update, param_vals, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60, log_params=None, which_organism='all', pro99_mode=False, t_eval=None):
    params = get_params(X, params_to_update, param_vals, log_params)
    return generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout, which_organism=which_organism, pro99_mode=pro99_mode, t_eval=t_eval)



def run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout=10*60, which_organism='all', pro99_mode=False, t_eval=None):
    try:
        run_args = [
             sys.executable, __file__, 
             '--json', json_fpath,
             '--ref_csv', ref_csv, 
             '--run_id', run_id,
             '--outdpath', out_dpath,
             '--which_organism', which_organism,

        ]
        if pro99_mode:
            run_args += ['--pro99_mode']
        if t_eval is not None:
            run_args += ['--t_eval'] + np.char.mod('%d', t_eval).tolist()

        #print(run_args)
        result = subprocess.run(run_args,
            capture_output=True, text=True, check=True, timeout=timeout,
        )
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)

        
    except subprocess.CalledProcessError as err:
        print('CalledProcessError', err.returncode)
        print("stdout:", err.stdout)
        print("stderr:", err.stderr)
                
    except subprocess.TimeoutExpired as err:
        print('TimeoutExpired', err.timeout)
        print("stdout:", err.stdout)
        print("stderr:", err.stderr)
    return run_id

def run_chunk(param_vals, param_values, params_to_update, chunk, number_of_runs, run_id, ref_csv, json_dpath, out_dpath, timeout, skip_if_found=True, log_params=None, which_organism='all', pro99_mode=False):
    start_line = (chunk  - 1) * number_of_runs
    end_line = min(param_values.shape[0], start_line + number_of_runs)
    if start_line >= end_line:
        print(f'start ({start_line}) >= end ({end_line}), stopping...')
        return
    print (f'running chunk {chunk} \t\t<<<==========')
    for i in range(start_line, end_line):

        out_fprefix = f'{run_id}_{i}'
        print(out_fprefix)

        files = glob.glob(os.path.join(out_dpath, f'{out_fprefix}_*_df.csv.gz'))
        if len(files) == 0:
            generate_json_and_run_from_X(
                param_values[i], params_to_update, param_vals, 
                ref_csv, json_dpath, out_dpath, out_fprefix, timeout, log_params=log_params, which_organism=which_organism, pro99_mode=pro99_mode)
            

def run_sensitivity_per_parameter(param_vals, parameter, bound, number_of_runs, run_id, ref_csv, json_dpath, out_dpath, timeout, skip_if_found=True, log_param=False, which_organism='all', pro99_mode=False):
    if log_param == True:
        bound = (np.log(bound[0]), np.log(bound[1]))
    for i,v in enumerate(np.linspace(bound[0], bound[1],num=number_of_runs)):
        out_fprefix = f'{run_id}_{parameter}_{i}'
        print(out_fprefix)
        generate_json_and_run_from_X(
            [v], [parameter], param_vals, 
            ref_csv, json_dpath, out_dpath, out_fprefix, timeout, log_params=[log_param], which_organism=which_organism, pro99_mode=pro99_mode)


    
if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', default='None')
    parser.add_argument('--json', help='json with param vals', nargs="+")
    parser.add_argument('--maxday', help='max day of simulation', type=int, default=140)

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)
    parser.add_argument("--model", help="model to run", choices=['MIN', 'FULL', 'LEAK', 'MIXO'], default='FULL')
    parser.add_argument("--which_organism", help="which organism to run", choices=['ponly', 'honly', 'all'], default='all')
    parser.add_argument("--pro99_mode", help="run on pro99 media",
                        action="store_true")
    parser.add_argument('--t_eval', nargs="+", type=float, default=None)
    
    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    if args.ref_csv == 'None':
        refdf = None
    else:
        refdf = pd.read_excel(args.ref_csv)
    param_vals = get_param_vals(args.model)

    #print(args.t_eval) 
    MSE_err = run_with_params_json(args.json, args.maxday, refdf, dpath, args.run_id, args.which_organism, args.pro99_mode, args.t_eval)
    print ('\nMSE:', MSE_err)

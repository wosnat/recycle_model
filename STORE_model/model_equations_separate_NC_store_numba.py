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

from numba import jit, njit



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
air_water_exchange_constant = (h / (Kg * B * 0.01))
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


DIN_IDX = 0
DON_IDX = 1
DIC_IDX = 2
DOC_IDX = 3
    
    
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
    param_df = pd.read_excel( 'Model_Parameters Store model.xlsx',)
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



def prepare_params_tuple_cc(param_vals):
    pars = param_vals
    paramCN  = np.array([pars['Rp'], pars['Rh']], dtype=np.float64)
    paramQCmax  = np.array([pars['QCmaxp'], pars['QCmaxh']], dtype=np.float64)
    paramQCmin  = np.array([pars['QCminp'], pars['QCminh']], dtype=np.float64)
    paramKmtb  = np.array([pars['Kmtbp'], pars['Kmtbh']], dtype=np.float64)
    paramb  = np.array([pars['bp'], pars['bh']], dtype=np.float64)
    paramr0  = np.array([pars['r0p'], pars['r0h']], dtype=np.float64)
    paramM = np.array([pars['Mp'], pars['Mh']], dtype=np.float64)
    paramOverflow = pars['OverflowMode']
    paramE_leak = np.array([pars['E_leakp'], pars['E_leakh']], dtype=np.float64)
    paramE_ROS = np.array([pars['E_ROSp'], pars['E_ROSp']], dtype=np.float64)
    paramomega_ROS = np.array([pars['omegaP'], pars['omegaH']], dtype=np.float64)
    
    paramVmaxROSh = pars['VmaxROSh']
    paramK_ROSh = pars['K_ROSh']
    paramROSMode = pars['ROSMode']
    paramgamma_DON2DIN = np.array([pars['gamma_DON2DINp'], pars['gamma_DON2DINh']], dtype=np.float64)
    paramgammaD = np.array([pars['gammaDp'], pars['gammaDh']], dtype=np.float64)
    paramROS_decay = pars['ROS_decay']
    kns = np.array(
        [[pars['KINp'], pars['KONp'], pars['KICp'], pars['KOCp']],
        [pars['KINh'], pars['KONh'], pars['KICh'], pars['KOCh'], ]], dtype=np.float64
    ).T
    vmax = np.array(
        [[pars['VmaxINp'], pars['VmaxONp'], pars['VmaxICp'], pars['VmaxOCp']],
        [pars['VmaxINh'], pars['VmaxONh'], pars['VmaxICh'], pars['VmaxOCh'], ]], dtype=np.float64
    ).T
    return (
        paramCN, paramQCmax, paramQCmin, kns, vmax, paramKmtb,paramOverflow, 
        paramb, paramr0, paramM, paramE_leak, paramE_ROS, paramVmaxROSh, paramK_ROSh,
        paramgamma_DON2DIN, paramgammaD, paramROS_decay, paramROSMode, paramomega_ROS)

def prepare_params_tuple_ponly(param_vals):
    pars = param_vals
    paramCN  = np.array([pars['Rp'], ], dtype=np.float64)
    paramQCmax  = np.array([pars['QCmaxp'], ], dtype=np.float64)
    paramQCmin  = np.array([pars['QCminp'], ], dtype=np.float64)
    paramKmtb  = np.array([pars['Kmtbp'], ], dtype=np.float64)
    paramb  = np.array([pars['bp'], ], dtype=np.float64)
    paramr0  = np.array([pars['r0p'], ], dtype=np.float64)
    paramM = np.array([pars['Mp'], ], dtype=np.float64)
    paramOverflow = pars['OverflowMode']
    paramE_leak = np.array([pars['E_leakp'], ], dtype=np.float64)
    paramE_ROS = np.array([pars['E_ROSp'], ], dtype=np.float64)
    paramomega_ROS = np.array([pars['omegaP'], ], dtype=np.float64)
    
    paramVmaxROSh = pars['VmaxROSh']
    paramK_ROSh = pars['K_ROSh']
    paramROSMode = pars['ROSMode']
    paramgamma_DON2DIN = np.array([pars['gamma_DON2DINp'], ], dtype=np.float64)
    paramgammaD = np.array([pars['gammaDp'], ], dtype=np.float64)
    paramROS_decay = pars['ROS_decay']
    kns = np.array(
        [[pars['KINp'], pars['KONp'], pars['KICp'], pars['KOCp']],
        ], dtype=np.float64
    ).T
    vmax = np.array(
        [[pars['VmaxINp'], pars['VmaxONp'], pars['VmaxICp'], pars['VmaxOCp']],
        ], dtype=np.float64
    ).T
    return (
        paramCN, paramQCmax, paramQCmin, kns, vmax, paramKmtb,paramOverflow, 
        paramb, paramr0, paramM, paramE_leak, paramE_ROS, paramVmaxROSh, paramK_ROSh,
        paramgamma_DON2DIN, paramgammaD, paramROS_decay, paramROSMode, paramomega_ROS)


@njit
def prepare_vals_tuple_cc(var_vals):
    Bp,   Np,    Cp,    Bh,    Nh,    Ch,     DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS = var_vals
    
    resources = np.array(
        [DIN, DON, DIC, DOC, ], dtype=np.float64
    ).reshape(4, 1)
    biomass = np.array([Bp, Bh], dtype=np.float64)
    storeN  = np.array([Np, Nh], dtype=np.float64)
    storeC  = np.array([Cp, Ch], dtype=np.float64)
    return (biomass, storeN, storeC, resources, DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS, Bh)

@njit
def prepare_vals_tuple_ponly(var_vals):
    Bp,   Np,    Cp,    DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS = var_vals
    Bh = 0 # for ROS
    resources = np.array(
        [DIN, DON, DIC, DOC, ], dtype=np.float64
    ).reshape(4, 1)
    biomass = np.array([Bp], dtype=np.float64)
    storeN  = np.array([Np], dtype=np.float64)
    storeC  = np.array([Cp], dtype=np.float64)
    return (biomass, storeN, storeC, resources, DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS, Bh)


@njit
def prepare_jac_sparsity_cc(ROS_mode):
    #element (i, j) is equal to d f_i / d y_j
    #Bp,Np,Cp,Bh,Nh,Ch,DON,RDON,DIN,DOC,RDOC,DIC,ROS 
    ROS = 1 if ROS_mode else 0
    return (np.array([
        #Bp,Np,Cp,Bh,Nh,Ch,DON,RDON,DIN,DOC,RDOC,DIC,ROS 
        [1 ,1, 1, 0, 0 ,0 ,1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Bp
        [1 ,1, 1, 0, 0 ,0 ,1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Np
        [1 ,1, 1, 0, 0 ,0 ,1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Cp
        [0, 0 ,0 ,1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Bh
        [0, 0 ,0 ,1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Nh
        [0, 0 ,0 ,1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # Ch
        [1 ,1, 1, 1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # DON
        [1 ,0 ,0 ,1 ,0, 0, 0  ,0   ,0  ,0  , 0  ,0  ,0  ],    # RDON
        [1 ,1, 1, 1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # DIN
        [1 ,1, 1, 1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # DOC
        [1 ,0 ,0 ,1 ,0, 0, 0  ,0   ,0  ,0  , 0  ,0  ,0  ],    # RDOC
        [1 ,1, 1, 1 ,1, 1, 1  ,0   ,1  ,1  , 0  ,1  ,ROS],    # DIC
        [1 ,0 ,0 ,1 ,0, 0, 0  ,0   ,0  ,0  , 0  ,0  ,1  ],    # ROS
        ]))
        
    


@njit
def compute_gross_uptake(
    biomass, storeN, storeC, resources, ROS,
    paramCN, paramQCmax, paramQCmin, paramkns, paramvmax, paramROSMode, paramomega_ROS, 
):
    QC = (storeC + biomass * paramCN) / (storeN + biomass)
  
    # monod ratios
    limits = resources / (resources + paramkns)
    # regulate uptake by not letting the stores grow too large 
    # question: is DIC uptake by photosynthesis regulated?
    regC = 1 - (QC / paramQCmax)
    regC = np.atleast_2d(regC)
    regN = 1 - (paramQCmin / QC)
    regN = np.atleast_2d(regN)
    reg =  np.vstack((regN, regN, regC, regC))
    reg_clipped = np.clip(reg, a_min =0.0, a_max=1.0)
    
    # ROS
    if paramROSMode == 1:
        ROS_penalty = np.exp( - paramomega_ROS*ROS)
    else:
        ROS_penalty = np.ones_like(paramomega_ROS)

    # gross uptake (regardless of C:N ratio)
    # vmax = muinfp* VmaxIp / Qp
    # umol N /L  or umol C /L  
    gross_uptake = paramvmax * biomass * limits  * reg_clipped * ROS_penalty
    uptakeN = gross_uptake[DIN_IDX,:] + gross_uptake[DON_IDX,: ]
    uptakeC = gross_uptake[DIC_IDX,:] + gross_uptake[DOC_IDX,: ]   
    return (gross_uptake, uptakeN, uptakeC, QC)

@njit
def compute_net_uptake(
    biomass, storeN, storeC, uptakeN, uptakeC,
    paramCN, paramKmtb, paramb, paramr0
):
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
    return(biosynthesisN, respirationC, biomass_breakdownC, netDeltaN, netDeltaC)

@njit
def compute_losses(
    biomass, storeN, storeC, netDeltaN, netDeltaC, paramCN, paramOverflow, paramM, paramE_leak
):
    # overflow -
    # make the store maintain the C:N ratio and exude the rest
    # Overflow quantity 
    # Oh/Op: enable overflow (0 or 1)
    # umol N / L
    if (paramOverflow == 1):
        store_keepN = np.minimum(netDeltaN, netDeltaC / paramCN) 
        overflowN = netDeltaN - store_keepN
        overflowC = netDeltaC - store_keepN * paramCN
    else:
        overflowN = np.zeros_like(netDeltaN)
        overflowC = np.zeros_like(netDeltaN)
        #overflowN = 0.0
        #overflowC = 0.0

    # death
    # Need to explain why we used exponential decay – in ISMEJ we show that other formulations are better for co-cultures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.
    deathN = paramM * biomass
    # leakiness formulated as fraction of biomass (“property tax”)
    # We assume that the vast majority of C and N biomass is in organic form, hence leakiness is to organic. We assume that overflow is also to organic in both organisms, as for the phototroph this is the release of fixed C (or inorganic N incorporated into e.g. AA) which cannot be used for growth. For the heterotrophs we assume overflow metabolism to be the inefficient use of organic C (e.g. not fully oxidized) to maximize growth rate (*citation E coli).
    leakinessN = paramE_leak * biomass
    return (deathN, leakinessN, overflowN, overflowC)

@njit
def compute_ROS(
    biomass, Bh, ROS, paramROSMode, paramROS_decay, paramE_ROS, paramVmaxROSh, paramK_ROSh
):
    if paramROSMode == 1:
        ROSdecay = ROS * paramROS_decay
        netROS = ROS - ROSdecay
        # ROS production depends on biomass
        ROSrelease = paramE_ROS * biomass
        ROSbreakdownh = paramVmaxROSh * Bh * netROS / (netROS + paramK_ROSh)
        dROSdt = np.sum(ROSrelease) - ROSdecay - ROSbreakdownh
        return dROSdt 
    else:
        return 0.0


@njit
def compute_C_odes(
    DIC, netDeltaC, gross_uptake, respirationC, overflowC, 
    deathN, leakinessN, paramCN, paramgammaD,  
):
    deathC = deathN * paramCN
    dic_air_water_exchange   = - (DIC - c_sat) / air_water_exchange_constant
    dCdt = netDeltaC - overflowC
    dDICdt = dic_air_water_exchange + np.sum(respirationC - gross_uptake[DIC_IDX,:])
    dDOCdt = np.sum(
        deathC * paramgammaD + leakinessN * paramCN + overflowC - gross_uptake[DOC_IDX,:])
    dRDOCdt = np.sum(deathC * (1 - paramgammaD))
    return(dCdt, dDICdt, dDOCdt, dRDOCdt)

@njit
def compute_N_odes(
    DON, biomass, netDeltaN, gross_uptake, biosynthesisN, biomass_breakdownC, overflowN, 
    deathN, leakinessN, paramCN, paramgammaD, paramgamma_DON2DIN
):
    # DON breakdown due to exoenzymes
    DON2DIN = paramgamma_DON2DIN * biomass * DON
    dBdt = biosynthesisN - biomass_breakdownC / paramCN - deathN - leakinessN 
    dNdt = netDeltaN - overflowN
    dDONdt = np.sum(
        deathN * paramgammaD + leakinessN - gross_uptake[DON_IDX,:] - DON2DIN
    )
    dRDONdt = np.sum(deathN * (1 - paramgammaD))

    # In discussion can state that if DIN is produced also through overflow or leakiness then this could support Pro growth, but this is not encoded into our model.
    # Assuming that recalcitrant DON isd released only during mortality (discuss release through leakiness)
    # Assuming RDON/RDOC is recalcitrant to both organisms   
    dDINdt = np.sum(overflowN + DON2DIN - gross_uptake[DIN_IDX,:])
    return(dBdt, dNdt, dDINdt, dDONdt, dRDONdt, DON2DIN)


#@njit
def basic_model_cc_ode_jit1(
    biomass, storeN, storeC, resources, DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS, Bh,
    paramCN, paramQCmax, paramQCmin, paramkns, paramvmax, paramKmtb,paramOverflow,
    paramb, paramr0, paramM, paramE_leak, paramE_ROS, paramVmaxROSh, paramK_ROSh,
    paramgamma_DON2DIN, paramgammaD, paramROS_decay, paramROSMode, paramomega_ROS
):


    gross_uptake, uptakeN, uptakeC, QC = compute_gross_uptake(
        biomass, storeN, storeC, resources, ROS,
        paramCN, paramQCmax, paramQCmin, paramkns, paramvmax, paramROSMode, paramomega_ROS)

    biosynthesisN, respirationC, biomass_breakdownC, netDeltaN, netDeltaC = compute_net_uptake(
        biomass, storeN, storeC, uptakeN, uptakeC,
        paramCN, paramKmtb, paramb, paramr0)

    deathN, leakinessN, overflowN, overflowC = compute_losses(
        biomass, storeN, storeC, netDeltaN, netDeltaC, paramCN, paramOverflow, paramM, paramE_leak)

    # final differential equations

    dCdt, dDICdt, dDOCdt, dRDOCdt = compute_C_odes(
        DIC, netDeltaC, gross_uptake, respirationC, overflowC, 
        deathN, leakinessN, paramCN, paramgammaD)

    dBdt, dNdt, dDINdt, dDONdt, dRDONdt, DON2DIN = compute_N_odes(
        DON, biomass, netDeltaN, gross_uptake, biosynthesisN, biomass_breakdownC, overflowN, 
        deathN, leakinessN, paramCN, paramgammaD, paramgamma_DON2DIN)
        
    dROSdt = compute_ROS(
        biomass, Bh, ROS, paramROSMode, paramROS_decay, paramE_ROS, paramVmaxROSh, paramK_ROSh)

    return (dBdt, dNdt, dCdt, 
        dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt,
        gross_uptake, uptakeN, uptakeC, QC,
        biosynthesisN, respirationC, biomass_breakdownC,
        overflowN, overflowC)


def basic_model_cc_ode(time, var_vals, par_tuple, return_intermediate=False):
    var_tuple = prepare_vals_tuple_cc(var_vals)
    (dBdt, dNdt, dCdt, 
        dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, 
        gross_uptake, uptakeN, uptakeC, QC,
        biosynthesisN, respirationC, biomass_breakdownC,
        overflowN, overflowC,
    ) = basic_model_cc_ode_jit1(*var_tuple,*par_tuple)
    if not return_intermediate:
        return (
            dBdt[0], dNdt[0], dCdt[0], 
            dBdt[1], dNdt[1], dCdt[1], 
            dDONdt, dRDONdt, dDINdt, 
            dDOCdt, dRDOCdt, dDICdt, dROSdt
        )
    else:
        intermediates = (
            gross_uptake, uptakeN, uptakeC, QC,
            biosynthesisN, respirationC, biomass_breakdownC,
            overflowN, overflowC,
        )
        flat_intermediates = [i  for f in intermediates for i in (f.flat if isinstance(f, np.ndarray) else [f])]

        return flat_intermediates


def basic_model_ponly_ode(time, var_vals, par_tuple, return_intermediate=False):
    var_tuple = prepare_vals_tuple_ponly(var_vals)
    (dBdt, dNdt, dCdt, 
        dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, 
        gross_uptake, uptakeN, uptakeC, QC,
        biosynthesisN, respirationC, biomass_breakdownC,
        overflowN, overflowC,
    ) = basic_model_cc_ode_jit1(*var_tuple,*par_tuple)
    if not return_intermediate:
        return (
            dBdt[0], dNdt[0], dCdt[0], 
            dDONdt, dRDONdt, dDINdt, 
            dDOCdt, dRDOCdt, dDICdt, dROSdt
        )
        
    else:
        intermediates = (
            gross_uptake, uptakeN, uptakeC, QC,
            biosynthesisN, respirationC, biomass_breakdownC,
            overflowN, overflowC,
        )
        flat_intermediates = [i  for f in intermediates for i in (f.flat if isinstance(f, np.ndarray) else [f])]

        return flat_intermediates
    

# functions

# number of cells (Q is a parameter)
# cells/L
#Xp = Bp / Qp
#Xh = Bh / Qh


def get_main_init_vars(pro99_mode):
    '''return : var_names, init_vars, intermediate_names'''
    var_names  = ['Bp', 'Np',  'Cp',  'Bh',  'Nh',  'Ch',   'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS',]
    intermediate_names = [
            'gross_uptakeINp', 'gross_uptakeINh', 'gross_uptakeONp', 'gross_uptakeONh', 
            'gross_uptakeICp', 'gross_uptakeICh', 'gross_uptakeOCp', 'gross_uptakeOCh', 
            'uptakeNp', 'uptakeNh', 'uptakeCp', 'uptakeCh', 'QCp','QCh',
            'biosynthesisNp', 'biosynthesisNh', 'respirationCp', 'respirationCh', 'biomass_breakdownCp','biomass_breakdownCh',
            'overflowNp', 'overflowNh', 'overflowCp','overflowCh',
    ]
    if pro99_mode:        
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_BH,INIT_NH,INIT_CH,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
    else:
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_BH,INIT_NH,INIT_CH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
    return var_names, np.array(init_vars, dtype=np.float64), intermediate_names


def get_ponly_init_vars(pro99_mode):
    '''return : var_names, init_vars, intermediate_names'''
    var_names  = ['Bp',  'Np',  'Cp',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    intermediate_names = [
            'gross_uptakeINp',  'gross_uptakeONp', 
            'gross_uptakeICp',  'gross_uptakeOCp', 
            'uptakeNp', 'uptakeCp', 'QCp',
            'biosynthesisNp',  'respirationCp',  'biomass_breakdownCp',
            'overflowNp', 'overflowCp',
    ]
    if pro99_mode:        
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS]
    else:
        init_vars = [INIT_BP,INIT_NP,INIT_CP,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS]
    return var_names, np.array(init_vars, dtype=np.float64), intermediate_names
    
def get_honly_init_vars(pro99_mode):
    '''return : var_names, init_vars, intermediate_names'''
    var_names  = ['Bh',  'DON', 'Nh',  'Ch',   'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    intermediate_names = [
            'gross_uptakeINh',  'gross_uptakeONh', 
            'gross_uptakeICh',  'gross_uptakeOCh', 
            'uptakeNh', 'uptakeCh', 'QCh',
            'biosynthesisNh',  'respirationCh',  'biomass_breakdownCh',
            'overflowNh', 'overflowCh',
    ]
    if pro99_mode:        
        init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN_PRO99,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
    else:
        init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
    return var_names, np.array(init_vars, dtype=np.float64), intermediate_names
    

def print_dydt0(calc_dydt, var_names, init_vars, par_tuple):
    dydt0 = calc_dydt(0, init_vars, par_tuple)
    for i,j, k in zip(var_names, dydt0, init_vars):
        print(f'd{i}/dt = {j:.2e}, init {i} = {k:.2e}, newval = {k+j:.2e}')


def print_intermediate0(intermediate_func, interm_names, init_vars):
    for i,j in zip(interm_names, intermediate_func(*init_vars)):
        print(f'{i:<4} = {j:.2e}')


def biomass_diff0(calc_dydt, var_names, init_vars, par_tuple):
    dydt0 = calc_dydt(0, init_vars, par_tuple)
    paramCN = par_tuple[0]
    V = dict(zip(var_names, dydt0))
    print (f"dBp/dt + dBh/dt + dDON/dt + dRDON/dt + dDIN/dt = { V['Bp'] + V['Np'] + V['Bh'] + V['Nh'] + V['DON'] + V['RDON'] + V['DIN'] }")
    print (f"dBp/dt + dBh/dt + dDOC/dt + dRDOC/dt + dDIC/dt = { V['Bp']*paramCN[0] + V['Cp'] + V['Bh']*paramCN[1] + V['Ch'] + V['DOC'] + V['RDOC'] + V['DIC'] }")

def biomass_diff0_ponly(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    V = dict(zip(var_names, dydt0))
    print (f"dBp/dt  + dDON/dt + dRDON/dt + dDIN/dt = { V['Bp']  + V['DON'] + V['RDON'] + V['DIN'] }")

def biomass_diff0_honly(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    V = dict(zip(var_names, dydt0))
    print (f"dBh/dt + dDON/dt + dRDON/dt + dDIN/dt = {  V['Bh'] + V['DON'] + V['RDON'] + V['DIN'] }")


def get_t_end(maxday=140, t_eval = None, seconds_in_day=seconds_in_day):
    
    if t_eval is None: 
        t_end = maxday*seconds_in_day
    else:
        t_end= np.max(t_eval)
    return t_end


def run_solver(calc_dydt, init_vars, par_tuple, t_end, t_eval, method='BDF', jac_sparsity=None):
    t_start = 0
    sol = solve_ivp(
        fun=calc_dydt, y0=init_vars, args=(par_tuple,),
        t_span=[t_start, t_end], t_eval=t_eval, max_step=1000, #first_step=1, 
        method=method, jac_sparsity=jac_sparsity)
        #method='Radau',)
    return sol

def solver2df(sol, var_names, par_tuple, t_eval=None,  intermediate_names=None, calc_dydt=None, ):
    d = dict(zip(var_names, sol.y))
    d['t'] = sol.t
    df = pd.DataFrame(data=d)
    if t_eval is not None:
        df = df.loc[np.rint(df.t).isin(np.rint(t_eval))].copy()
    df['day'] = df['t']/seconds_in_day
    if intermediate_names: 
        df[intermediate_names] = df[var_names].apply(
            lambda x : calc_dydt(0, x.to_numpy(), par_tuple=par_tuple, return_intermediate=True), axis=1, 
            result_type='expand')
    paramCN = par_tuple[0]
    if 'Bp' in df.columns:
        df['Bp[C]'] = df['Bp']*paramCN[0]
        df['Bptotal'] = df['Bp']+df['Np']
        df['Bptotal[C]'] = df['Bp[C]']+df['Cp']
        
    if 'Bh' in df.columns:
        # TODO: will not work in honly mode
        df['Bh[C]'] = df['Bh']*paramCN[1]
        df['Bhtotal'] = df['Bh']+df['Nh']
        df['Bhtotal[C]'] = df['Bh[C]']+df['Ch']
        
    return df

   
def _mse(x, refdf, refcol, col, timecol, tolerance):
    #print(x.columns)
    tdf = pd.merge_asof(x[[timecol, col]].dropna(), refdf[[timecol, refcol]], on=timecol, direction='nearest', tolerance=tolerance).dropna()
    return pd.Series({
        'compare_points': tdf.shape[0], 
        'MSE': mean_squared_error(tdf[col], tdf[refcol])}
    )

def compute_mse(df, refdf, refcol= 'ref_Bp', col='Bptotal', timecol='t', tolerance=100):
    #print (x)
    mse_df = refdf.groupby(['Sample', 'full name', 'Group',]
                        ).apply(lambda y : _mse(df,y, refcol= refcol, col=col, timecol=timecol, tolerance=tolerance))
    mse_df = mse_df.reset_index()
    return mse_df

def run_solver_from_new_params(
    new_param_vals, refdf, init_var_vals, 
    calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names, return_dfs=False
    ):

    par_tuple = prepare_params_tuple(new_param_vals)
    sol = run_solver(calc_dydt, init_var_vals, par_tuple, t_end , t_eval)
    df = solver2df(
        sol, var_names, par_tuple=par_tuple, 
        intermediate_names=intermediate_names, calc_dydt=calc_dydt)

    if refdf is not None:
        try:
            mse_df = compute_mse(df, refdf, refcol= 'ref_Bp', col='Bptotal', timecol='t', tolerance=100)
            perr = mse_df['MSE'].min()

        except Exception as inst:
            print(inst)
            pass
    if return_dfs:
        return(perr, new_param_vals, df, mse_df)
    else:
        return perr 

def save_solver_results_to_file(
    new_param_vals, df, mse_df,
    out_dpath, out_fprefix, 
    ):
        
    sumdf = pd.DataFrame({str(k): v for k,v in new_param_vals.items()}, index=[0])
    sumdf['run_id'] = out_fprefix
    df['run_id'] = out_fprefix
    mse_df['run_id'] = out_fprefix
    
    df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'), compression='gzip', index=False)
    sumdf.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'), compression='gzip', index=False)
    mse_df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_mse.csv.gz'), compression='gzip', index=False)



def run_solver_from_X(
    X, params_to_update, orig_param_vals, refdf, 
    log_params, init_var_vals, 
    calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names, return_dfs=False
    ):

    new_param_vals = get_params(X, params_to_update, orig_param_vals, log_params)
    return run_solver_from_new_params(
        new_param_vals, refdf, init_var_vals, 
        calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names, return_dfs
    )


def run_solver_from_new_params_and_save(
    new_param_vals, refdf, out_dpath, 
    out_fprefix, init_var_vals, 
    calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names,
    ):

    perr, new_param_vals, df, mse_df,  =  run_solver_from_new_params(
        new_param_vals, refdf, init_var_vals, 
        calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names, return_dfs=True
    )
            
    save_solver_results_to_file(
        new_param_vals, df, mse_df,
        out_dpath, out_fprefix, 
    )

    return perr 



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

def get_params(X, params_to_update, param_vals, log_params=None):
    new_param_vals = param_vals.copy()
    if log_params is not None:
        X = [np.exp(x) if lg else x for x,lg in zip(X, log_params)]

    new_param_vals.update({k : v for k,v in zip(params_to_update, X)})
    return new_param_vals




def get_runid_unique_suffix(pro99_mode, param_vals):
    suffix = ''
    if pro99_mode:
        suffix = '_pro99'
    hash_val = str(hash(tuple(param_vals.values())))
    suffix = f'{suffix}_h{hash_val}'
    return suffix 

    
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
    parser.add_argument("--model", help="model to run", choices=['MIN', 'MIXOTROPHY', 'OVERFLOW', 'ROS', 'EXOENZYMES'], default='FULL')
    parser.add_argument("--which_organism", help="which organism to run", choices=['ponly', 'honly', 'all'], default='all')
    parser.add_argument("--pro99_mode", help="run on pro99 media",
                        action="store_true")
    parser.add_argument('--t_eval', nargs="+", type=float, default=None)
    parser.add_argument("--param_sensitivity", help="index of param to update (0 based) ",  type=int, default=-1)
    parser.add_argument("--organism_to_tune", help="which organism to tune", choices=['PRO', 'HET'], default='PRO')
    parser.add_argument("--number_of_runs", help="number of simulations to run",  type=int, default=1024)
    
    
    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)

    if args.ref_csv == 'None':
        refdf = None
    else:
        refdf = pd.read_excel(args.ref_csv)

    new_param_vals = get_param_vals(args.model)
    if args.json is not None:
        for json_fpath in args.json:
            new_param_vals = json2params(new_param_vals, json_fpath)

    suffix = get_runid_unique_suffix(args.pro99_mode, new_param_vals)
        
    t_eval = args.t_eval
    if t_eval is None:
        if refdf is not None:
            t_eval = np.rint(refdf['t'].drop_duplicates().sort_values()).to_numpy()
    t_end = get_t_end(maxday=args.maxday, t_eval=t_eval, seconds_in_day=seconds_in_day)

    if t_eval is None: # no refdf
        t_eval = np.arange(0, t_end, 1000)

    which_organism = args.which_organism
    if which_organism == 'ponly':
        var_names, init_var_vals, intermediate_names =  get_ponly_init_vars(pro99_mode=args.pro99_mode)
        calc_dydt = basic_model_ponly_ode
        prepare_params_tuple = prepare_params_tuple_ponly
    elif which_organism == 'honly':
        var_names, init_var_vals, intermediate_names =  get_honly_init_vars(pro99_mode=args.pro99_mode)
        # TODO - HONLY implementation
        calc_dydt = basic_model_ponly_ode
        prepare_params_tuple = prepare_params_tuple_ponly
    else:
        var_names, init_var_vals, intermediate_names =  get_main_init_vars(pro99_mode=args.pro99_mode)
        calc_dydt = basic_model_cc_ode
        prepare_params_tuple = prepare_params_tuple_cc

    if args.param_sensitivity != -1:
        # run sensitivity
        idx = args.param_sensitivity
        params_to_update, bounds, log_params = get_param_tuning_values(args.model, args.organism_to_tune)
        parameter = params_to_update[idx]
        param_bounds = bounds[idx]
        param_log = log_params[idx]
        if param_log:
            param_bounds = (np.log(param_bounds[0]), np.log(param_bounds[1]))
        number_of_runs = args.number_of_runs
        for i,v in enumerate(np.linspace(param_bounds[0], param_bounds[1],num=number_of_runs)):
            if param_log:
                v = np.exp(v)
            new_param_vals[parameter] = v
            run_id = f"{args.run_id}_{args.model}_{parameter}_{i}_{v}_{suffix}"
            print(run_id)
            
            MSE_err = run_solver_from_new_params_and_save(
                new_param_vals, refdf, args.outdpath, 
                run_id, init_var_vals, 
                calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names,
            )
            print ('MSE:', MSE_err)
            
    else:
        # default - run simulation
        run_id = f'{args.run_id}_{suffix}'
        MSE_err = run_solver_from_new_params_and_save(
            new_param_vals, refdf, args.outdpath, 
            run_id, init_var_vals, 
            calc_dydt, prepare_params_tuple, t_end , t_eval, var_names, intermediate_names,
        )

        print ('MSE:', MSE_err)
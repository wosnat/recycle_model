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

# variables
Bp, Bh, DOC, RDOC, DIC, DON, RDON, DIN, ROS, ABp, ABh = symbols('B_p B_h DOC RDOC DIC DON RDON DIN ROS AB_p AB_h')


# parameters
gammaDp, gammaDh, EOp, EIp, EOh, EIh = symbols('gamma^D_p gamma^D_h E^O_p E^I_p E^O_h E^I_h')
KABp, KABh, EABp, EABh, MABp, MABh, decayABp, decayABh = symbols('K^AB_p K^AB_h E^AB_p E^AB_h M^AB_p M^AB_h decay^AB_p decay^AB_h')

Op, Oh = symbols('O_p O_h')
E_ROSp, E_ROSh, VROSmax, K_ROSh, omegaP, omegaH, ROS_decay= symbols('E^ROS_p E^ROS_h  VROSmax KROS_h omega_p omega_h ROS_decay')
Mp, Mh = symbols('M_p M_h')
Rp, Rh = symbols('R_p R_h')

KONp, KINp, KOCp, KICp, KONh, KINh, KOCh, KICh = symbols('K^ON_p K^IN_p K^OC_p K^IC_p K^ON_h K^IN_h K^OC_h K^IC_h')
VmaxONp, VmaxINp, VmaxOCp, VmaxICp, VmaxONh, VmaxINh, VmaxOCh, VmaxICh = symbols('Vmax^ON_p Vmax^IN_p Vmax^OC_p Vmax^IC_p Vmax^ON_h Vmax^IN_h Vmax^OC_h Vmax^IC_h')

tau, r0p, r0h, bp, bh = symbols('tau r0_p r0_h b_p b_h')



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


# initial concentrations
INIT_DIN = 100
INIT_DON = 20
INIT_RDON = 0
INIT_RDOC = 0
# Dalit DIC: 1618.825333  or 1.62E+03 uM
INIT_DIC = 3000

# Dalit: DOC 0.047155376888899 nM ?
# Dalit init TOC 16 mM
INIT_DOC = INIT_DON * R_CN
INIT_BP = 1e9 * Qp
INIT_BH = 1e10 * Qh
INIT_BH_CC = 5e9 * Qh # the actual concentration in the measurements
INIT_ROS = 0.2 # Morris, J. Jeffrey, et al. "Dependence of the cyanobacterium Prochlorococcus on hydrogen peroxide scavenging microbes for growth at the ocean's surface." PloS one 6.2 (2011): e16805.‏
INIT_SP = 0
INIT_SH = 0






# DIC exchange
h = 0.3 # height in m
Kg = 0.4968 / seconds_in_day # m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
c_sat = INIT_DIC


param_vals_with_symbols = {
    # Mortality and refractoriness
    # 1/d
    Mh : 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.6,         # pro death release 
    gammaDh : 0.6,         # het death release
    
    # C:N ratio
    Rp : R_P,
    Rh : R_H,
    
    # leakiness rates (“property tax”)
    # 1/d
    # EOp : 0.2 / seconds_in_day,        # move to income tax # pro organic leakiness 
    # EIh : 0.2 / seconds_in_day,        # move to income tax # het inorganic leakiness 
    EOp : 0.1 / seconds_in_day,        # pro organic leakiness 
    EIp : 0 / seconds_in_day,          # pro inorganic leakiness
    EOh : 0.1 / seconds_in_day,          # het organic leakiness 
    EIh : 0 / seconds_in_day,        # het inorganic leakiness 
    
    # K’s (affinity).
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake

    KONp : 0.17 * pro_vol**0.27, 
    KINp : 0.17 * pro_vol**0.27, # x 5 based on sensitivity
    KOCp : 0.17 * pro_vol**0.27, 
    KICp : 0.17 * pro_vol**0.27, 
    KONh : 0.17 * alt_vol**0.27, 
    KINh : 0.17 * alt_vol**0.27, # / 10 based on sensitivity
    KOCh : 0.17 * alt_vol**0.27, # x 5 based on sensitivity
    KICh : 0.17 * alt_vol**0.27, 

    # 1/d
    VmaxONp : 0 / seconds_in_day, 
    VmaxINp : 0.7 / seconds_in_day, 
    VmaxOCp : 0 * R_P /  seconds_in_day, 
    VmaxICp : 0.7 * R_P / seconds_in_day, 
    VmaxONh : 1.3 / seconds_in_day, 
    VmaxINh : 1.3 / seconds_in_day, # x 3 based on sensitivity
    VmaxOCh : 1.3 * R_H / seconds_in_day, 
    VmaxICh : 0 * R_H  / seconds_in_day, 
    
    #Vmaxp : 0.8 * 1.9e-9 * pro_vol**0.67 /Qp / seconds_in_day, 
    #Vmaxh : 3 * 1.9e-9 * alt_vol**0.67 / Qh / seconds_in_day, 

    
    # ROS

    # umol/cell/d
    E_ROSp : 1e-10 / Qp / seconds_in_day, 
    E_ROSh : 1e-10 / Qh / seconds_in_day, 
    VROSmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    K_ROSh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omegaP : 0.01, #0.1,
    omegaH : 0.00001, #0.1,
    
    # ROS decay 
    ROS_decay : 0.01,
    
    # antibiotic
    # umolN/L
    KABp : 0.17 * pro_vol**0.27 * 100, 
    KABh : 0.17 * pro_vol**0.27 * 100, 
    # umol/cell/d
    EABp : 0, #1e-20 / Qp / seconds_in_day, 
    EABh : 0, #1e-20 / Qh / seconds_in_day, 
    # antibiotic decay (1/d)
    decayABh : 0.01, 
    decayABp : 0.01, 
    
    # 1/d (between -1 to 1)  / seconds_in_day
    MABp : -0.01 / seconds_in_day, 
    MABh : -0.01 / seconds_in_day,
    # DIC (CO2) 
    tau : h / Kg,
    # dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989 
    r0p : 0.18 / seconds_in_day, 
    r0h : 0.18 / seconds_in_day,  
    # respiration coefficient, no units, Geider & Osborne 1989
    bp  : 0.01, 
    bh  : 0.01,
    
}
param_vals = {str(k) : v for k,v in param_vals_with_symbols.items()}


param_vals_neutral_with_symbols = {
    # 1/d
    Mh: 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.8,         # pro death release 
    gammaDh : 0.8,         # het death release
    
    # ratio
    Rp : R_P,
    Rh : R_H,
    
    # 1/d
    # EOp : 0.2 / seconds_in_day,        # move to income tax # pro organic leakiness 
    # EIh : 0.2 / seconds_in_day,        # move to income tax # het inorganic leakiness 
    EOp : 0.05 / seconds_in_day,        # pro organic leakiness 
    EIp : 0 / seconds_in_day,          # pro inorganic leakiness
    EOh : 0.05 / seconds_in_day,          # het organic leakiness 
    EIh : 0 / seconds_in_day,        # het inorganic leakiness 
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake

    KONp : 0.17 * pro_vol**0.27, 
    KINp : 0.17 * pro_vol**0.27, # x 5 based on sensitivity
    KOCp : 0.17 * pro_vol**0.27, 
    KICp : 0.17 * pro_vol**0.27, 
    KONh : 0.17 * alt_vol**0.27, 
    KINh : 0.17 * alt_vol**0.27, # / 10 based on sensitivity
    KOCh : 0.17 * alt_vol**0.27, # x 5 based on sensitivity
    KICh : 0.17 * alt_vol**0.27, 
    # umol N/cell/d
    # vmax = muinfp* VmaxIp * Qp
    # 1/day  * umol/cell  * umol/cell/d # TODO - figure out units
    VmaxONp : 0.7 * 1.9e-9 / Qp / seconds_in_day, 
    VmaxINp : 0.7 * 1.9e-9 / Qp / seconds_in_day, 
    VmaxOCp : 0.7 * 1.9e-9 * R_P / Qp / seconds_in_day, 
    VmaxICp : 0.7 * 1.9e-9 * R_P / Qp / seconds_in_day, 
    VmaxONh : 2 * 1.9e-9 / Qh / seconds_in_day, 
    VmaxINh : 2 * 1.9e-9 / Qh / seconds_in_day, # x 3 based on sensitivity
    VmaxOCh : 2 * 1.9e-9 * R_H / Qh / seconds_in_day, 
    VmaxICh : 2 * 1.9e-9 / 10000 * R_H / Qh / seconds_in_day, 
    
    #Vmaxp : 0.8 * 1.9e-9 * pro_vol**0.67 /Qp / seconds_in_day, 
    #Vmaxh : 3 * 1.9e-9 * alt_vol**0.67 / Qh / seconds_in_day, 

    
    # umol/cell/d
    E_ROSp : 0, # 1e-10 / seconds_in_day, 
    E_ROSh : 0, # 1e-10 / seconds_in_day, 
    VROSmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    K_ROSh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l 
    # How much do the ROS affect growth

    omegaP : 0.01, #0.1,
    omegaH : 0.00001, #0.1,
    
    ROS_decay : 0.01, 
    
    # Signals
    # umolN/L
    KABp : 0.17 * pro_vol**0.27 * 100, 
    KABh : 0.17 * pro_vol**0.27 * 100, 
    # umol/cell/d
    EABp : 0,  #1e-20 / Qp / seconds_in_day, 
    EABh : 0, # 1e-20 / Qh / seconds_in_day, 
    # 1/d (between -1 to 1)  / seconds_in_day
    MABp : 0.01 / seconds_in_day, 
    MABh : 0.01 / seconds_in_day,
    # antibiotic decay (1/d)
    decayABh : 0.01, 
    decayABp : 0.01, 

    # DIC (CO2) 
    tau : h / Kg,
    # dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989 
    r0p : 0.18 / seconds_in_day, 
    r0h : 0.18 / seconds_in_day,  
    # respiration coefficient, no units, Geider & Osborne 1989
    bp  : 0.01, 
    bh  : 0.01,
    
}
param_vals_neutral = {str(k) : v for k,v in param_vals_neutral_with_symbols.items()}



from enum import Flag, auto
class DISABLE_MECHANISMS(Flag):
    COMPETITION = auto()
    P_RECYCLING = auto()
    H_RECYCLING = auto()
    P_EXUDATION = auto()
    H_EXUDATION = auto()
    P_OVERFLOW = auto()
    H_OVERFLOW = auto()
    MIXOTROPHY = auto()
    DETOXIFICATION = auto()
    P_SIGNAL = auto()
    H_SIGNAL = auto()


def disable_mechanism(mechanisms, param_vals):
    ''' change the param vals to disable a given interaction mechanism
    use | to disable multiple mechanisms
    
    '''
    new_param_vals = param_vals.copy()
    if DISABLE_MECHANISMS.COMPETITION in mechanisms:
        new_param_vals[str(VmaxINh)] = new_param_vals[str(VmaxINh)] * 1e-3
    if DISABLE_MECHANISMS.P_RECYCLING in mechanisms:
        new_param_vals[str(gammaDp)] = 0
    if DISABLE_MECHANISMS.H_RECYCLING in mechanisms:
        new_param_vals[str(gammaDh)] = 0
    if DISABLE_MECHANISMS.P_EXUDATION in mechanisms:
        new_param_vals[str(EOp)] = 0
        new_param_vals[str(EIp)] = 0
    if DISABLE_MECHANISMS.H_EXUDATION in mechanisms:
        new_param_vals[str(EOh)] = 0
        new_param_vals[str(EIh)] = 0
    if DISABLE_MECHANISMS.P_OVERFLOW in mechanisms:
        new_param_vals[str(Op)] = 0
    if DISABLE_MECHANISMS.H_OVERFLOW in mechanisms:
        new_param_vals[str(Oh)] = 0
    if DISABLE_MECHANISMS.MIXOTROPHY in mechanisms:
        new_param_vals[str(VmaxONp)] = new_param_vals[str(VmaxONp)] * 50
        new_param_vals[str(VmaxOCp)] = new_param_vals[str(VmaxOCp)] * 50
    if DISABLE_MECHANISMS.DETOXIFICATION in mechanisms:
        new_param_vals[str(omegaP)] = 0
        new_param_vals[str(omegaH)] = 0
    if DISABLE_MECHANISMS.P_SIGNAL in mechanisms:
        new_param_vals[str(EABp)] = 0
    if DISABLE_MECHANISMS.H_SIGNAL in mechanisms:
        new_param_vals[str(EABh)] = 0

    return new_param_vals



# different model configurations
DISABLE_ROS_SIGNAL = {
    # ROS
    str(E_ROSh), 
    str(E_ROSp), 
    str(VROSmax), 
    str(K_ROSh), 
    str(omegaP), #0.1,
    str(omegaH), #0.1,
    
    # Signals
    # umolN/L
    str(KABp), 
    str(KABh), 
    # umol/cell/d
    str(EABp), 
    str(EABh), 
    # 1/d (between -1 to 1)  / seconds_in_day
    str(MABp), 
    str(MABh),
}

DISABLE_MIXOTROPHY = {
    str(VmaxONp), 
    str(VmaxOCp), 
}

DISABLE_EXUDATION_OVERFLOW = {
    str(EOp),
    str(EIp),
    str(EOh),
    str(EIh),
}


def model_partial(disable_set, param_vals, params_to_update, bounds, log_params):
    new_param_vals = param_vals.copy()
    # zero out the disabled params
    new_param_vals.update({str(p): 0 for p in disable_set})
    new_params_to_update = [p for p in params_to_update if p not in disable_set]
    new_bounds = [b for p,b in zip(params_to_update, bounds) if p not in disable_set]
    new_log_params = [b for p,b in zip(params_to_update, log_params) if p not in disable_set]
    return new_param_vals, new_params_to_update, new_bounds, new_log_params
    
def model_min(param_vals, params_to_update, bounds, log_params):
    return model_partial(
        DISABLE_ROS_SIGNAL | DISABLE_MIXOTROPHY | DISABLE_EXUDATION_OVERFLOW, 
        param_vals, params_to_update, bounds, log_params
    )

def model_mixotrophy(param_vals, params_to_update, bounds, log_params):
    return model_partial(
        DISABLE_ROS_SIGNAL, 
        param_vals, params_to_update, bounds, log_params
    )


def model_exudation(param_vals, params_to_update, bounds, log_params):
    return model_partial(
        DISABLE_ROS_SIGNAL | DISABLE_MIXOTROPHY, 
        param_vals, params_to_update, bounds, log_params
    )

def set_model(modeltype, param_vals, params_to_update, bounds, log_params):
    if modeltype == 'full':
        return param_vals, params_to_update, bounds, log_params
    elif modeltype == 'min':
        disableset = DISABLE_ROS_SIGNAL | DISABLE_MIXOTROPHY | DISABLE_EXUDATION_OVERFLOW
    elif modeltype == 'exu':
        disableset = DISABLE_ROS_SIGNAL | DISABLE_MIXOTROPHY 
    elif modeltype == 'mix':
        disableset = DISABLE_ROS_SIGNAL  
    return model_partial(disableset, 
        param_vals, params_to_update, bounds, log_params
    )
        
     

def print_params(param_vals=param_vals):
    for i in param_vals:
        print(i, f' = {param_vals[i]:.2e}, {param_vals[i] * seconds_in_day:.2e}')

# functions

# number of cells (Q is a parameter)
# cells/L
Xp = Bp / Qp
Xh = Bh / Qh


# monod ratios
limINp = (DIN / (DIN + KINp))
limONp = (DON / (DON + KONp))
limICp = (DIC / (DIC + KICp))
limOCp = (DOC / (DOC + KOCp))
limINh = (DIN / (DIN + KINh))
limONh = (DON / (DON + KONh))
limICh = (DIC / (DIC + KICh))
limOCh = (DOC / (DOC + KOCh))

limABh = (ABh / (ABh + KABh))
limABp = (ABp / (ABp + KABp))


# gross uptake (regardless of C:N ratio)
# vmax = muinfp* VmaxIp / Qp
# umol N /L 
gross_uptakeINp = VmaxINp * limINp * exp(-omegaP*ROS) * Bp
gross_uptakeONp = VmaxONp * limONp * exp(-omegaP*ROS) * Bp
gross_uptakeINh = VmaxINh * limINh * exp(-omegaH*ROS) * Bh
gross_uptakeONh = VmaxONh * limONh * exp(-omegaH*ROS) * Bh
# umol C /L
gross_uptakeICp = VmaxICp * limICp * exp(-omegaP*ROS) * Bp 
gross_uptakeOCp = VmaxOCp * limOCp * exp(-omegaP*ROS) * Bp 
gross_uptakeICh = VmaxICh * limICh * exp(-omegaH*ROS) * Bh 
gross_uptakeOCh = VmaxOCh * limOCh * exp(-omegaH*ROS) * Bh 

# net uptake (maintains C:N ratios)
# umol N / L
net_uptakeNp = Min(gross_uptakeINp + gross_uptakeONp, 
                     (gross_uptakeICp + gross_uptakeOCp) / Rp)

net_uptakeNh = Min(gross_uptakeINh + gross_uptakeONh, 
                     (gross_uptakeICh + gross_uptakeOCh) / Rh)

# inorganic to organic uptake ratio
# IOuptakeRateNp = gross_uptakeINp / (gross_uptakeINp + gross_uptakeONp)
# IOuptakeRateCp = gross_uptakeICp / (gross_uptakeICp + gross_uptakeOCp)
# IOuptakeRateNh = gross_uptakeINh / (gross_uptakeINh + gross_uptakeONh)
# IOuptakeRateCh = gross_uptakeICh / (gross_uptakeICh + gross_uptakeOCh)


# Overflow quantity 
# umol N / L
overflowNp = gross_uptakeINp + gross_uptakeONp - net_uptakeNp 
overflowNh = gross_uptakeINh + gross_uptakeONh - net_uptakeNh
# umol C /L
overflowCp = gross_uptakeICp + gross_uptakeOCp - net_uptakeNp * Rp
overflowCh = gross_uptakeICh + gross_uptakeOCh - net_uptakeNh * Rh

# umol N / L

# Respiration – growth associated bp/bh and maintenance associated r0p/r0h
# b * growth + r0 * biomass
respirationp =  bp* net_uptakeNp + Bp * r0p
respirationh =  bh* net_uptakeNh + Bh * r0h
dic_air_water_exchange   =  - (DIC - c_sat) / tau




# M = M / Q

# death = M * X = M * B /Q = M / Q * B
# death rate should be between 0 - 1
# Maximum death rate everyone dead in one sec
# minimum death rate - 0
# If only AB then effect of AB on mortality is only positive
death_ratep = Min(Mp + MABp*limABh, 1 )
death_rateh = Min(Mh + MABh*limABp, 1 )


# Need to explain why we used exponential decay – in ISMEJ we show that other formulations are better for co-cuiltures but these are emergent properties which we are explicitly testing here, and for the axenic cultures the exponential decay was good.

deathp = death_ratep * Bp 
deathh = death_rateh * Bh 

# leakiness formulated as fraction of biomass (“property tax”)
leakinessOp = EOp * Bp
leakinessIp = EIp * Bp
leakinessOh = EOh * Bh
leakinessIh = EIh * Bh

# ROS production depends on biomass
# epsilon = epsilon / Q
# VTMax = VROSmax / Q

ROSreleasep = E_ROSp * Bp
ROSreleaseh = E_ROSh * Bh
ROSbreakdownh = VROSmax * ROS / (ROS + K_ROSh) * Bh

# AB
# TODO Change to antibiotic
ABreleasep = EABp * Bp 
ABreleaseh = EABh * Bh

# We assume that the vast m,ajority of C and N biomass is in organic form, hence leakiness is to organic. We assume that overflow is also to organic in both organisms, as for the phototroph this is the release of fixed C (or inorganic N incorporated into e.g. AA) which cannot be used for growth. For the heterotrophs we assume overflow metabolism to be the inefficient use of organic C (e.g. not fully oxidized) to maximize growth rate (*citation E coli).

####################################################
# final equation - coculture
dBpdt = net_uptakeNp - deathp - leakinessOp - respirationp - ABreleasep
dBhdt = net_uptakeNh - deathh - leakinessOh - respirationh - ABreleaseh

# todo disable overflow
dDONdt = (
    deathp * gammaDp + leakinessOp + overflowNp - gross_uptakeONp +
    deathh * gammaDh + leakinessOh + overflowNh - gross_uptakeONh
    ) 
dDOCdt = (
    deathp * gammaDp * Rp + leakinessOp * Rp + overflowCp - gross_uptakeOCp +
    deathh * gammaDh * Rh + leakinessOh * Rh + overflowCh - gross_uptakeOCh)


# In discussion can state that if DIN is produced also through overflow or leakiness then this could support Pro growth, but this is not encoded into our model.
# Assuming that recalcitrant DON isd released only during mortality (discuss release through leakiness)
# Assuming RDON/RDOC is recalcitrant to both organisms

dRDONdt = (
    deathp * (1 - gammaDp) + 
    deathh * (1 - gammaDh))
dRDOCdt = (
    deathp * (1 - gammaDp) * Rp + 
    deathh * (1 - gammaDh) * Rh)

# Respiration of N is not a biological reality in this case (no NO3 respiration), and is used to maintain C:N ratio. It can be thought of as the release of NH4/urea for example during AA degradation

# Point for discussion with Mick 
dDINdt = (
    respirationp - gross_uptakeINp +
    respirationh - gross_uptakeINh)
# leakiness of inorganic
dDICdt = (
    dic_air_water_exchange +
    respirationp * Rp - gross_uptakeICp +
    respirationh * Rh - gross_uptakeICh)

# TODO Need to explain why Max and not just ROS dynamics
dROSdt = Max( -ROS*ROS_decay + ROSreleasep + ROSreleaseh - ROSbreakdownh, -ROS)
dABpdt = ABreleasep - ABp*decayABp
dABhdt = ABreleaseh - ABh*decayABh



####################################################
# PRO only model
dDONdt_ponly = (
    deathp * gammaDp + leakinessOp + overflowNp - gross_uptakeONp) 
dDOCdt_ponly = (
    deathp * gammaDp * Rp + leakinessOp * Rp + overflowCp - gross_uptakeOCp)


# In discussion can state that if DIN is produced also through overflow or leakiness then this could support Pro growth, but this is not encoded into our model.
# Assuming that recalcitrant DON isd released only during mortality (discuss release through leakiness)
# Assuming RDON/RDOC is recalcitrant to both organisms

dRDONdt_ponly = (
    deathp * (1 - gammaDp))
dRDOCdt_ponly = (
    deathp * (1 - gammaDp) * Rp)

# Respiration of N is not a biological reality in this case (no NO3 respiration), and is used to maintain C:N ratio. It can be thought of as the release of NH4/urea for example during AA degradation

# Point for discussion with Mick 
dDINdt_ponly = (
    respirationp - gross_uptakeINp)
# leakiness of inorganic
dDICdt_ponly = (
    dic_air_water_exchange +
    respirationp * Rp - gross_uptakeICp)

# TODO Need to explain why Max and not just ROS dynamics
dROSdt_ponly = Max( -ROS*ROS_decay + ROSreleasep, -ROS)
dABhdt_ponly = Integer(0)

####################################################
# HET only model
dDONdt_honly = (
    deathh * gammaDh + leakinessOh + overflowNh - gross_uptakeONh
    ) 
dDOCdt_honly = (
    deathh * gammaDh * Rh + leakinessOh * Rh + overflowCh - gross_uptakeOCh)


# In discussion can state that if DIN is produced also through overflow or leakiness then this could support Pro growth, but this is not encoded into our model.
# Assuming that recalcitrant DON isd released only during mortality (discuss release through leakiness)
# Assuming RDON/RDOC is recalcitrant to both organisms

dRDONdt_honly = (
    deathh * (1 - gammaDh))
dRDOCdt_honly = (
    deathh * (1 - gammaDh) * Rh)

# Respiration of N is not a biological reality in this case (no NO3 respiration), and is used to maintain C:N ratio. It can be thought of as the release of NH4/urea for example during AA degradation

# Point for discussion with Mick 
dDINdt_honly = (
    respirationh - gross_uptakeINh)
# leakiness of inorganic
dDICdt_honly = (
    dic_air_water_exchange +
    respirationh * Rh - gross_uptakeICh)

# TODO Need to explain why Max and not just ROS dynamics
dROSdt_honly = Max( -ROS*ROS_decay + ROSreleaseh - ROSbreakdownh, -ROS)
dABpdt_honly = Integer(0)


def print_equations():
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'ABp', 'ABh']
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dABpdt, dABhdt]
    for n,f in zip(var_names, sfunc_list):
        print(f'd{n}/dt')
        display(f)

def get_main_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dABpdt, dABhdt]
    var_list   = [ Bp,    Bh,    DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS,    ABp,    ABh]
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS',   'ABp',  'ABh']
    init_vars = [INIT_BP,INIT_BH_CC,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xp, Xh,

        limINp,
        limONp,
        limICp,
        limOCp,
        limINh,
        limONh,
        limICh,
        limOCh,

        gross_uptakeINp,
        gross_uptakeONp,
        gross_uptakeICp,
        gross_uptakeOCp,
        gross_uptakeINh,
        gross_uptakeONh,
        gross_uptakeICh,
        gross_uptakeOCh,

        net_uptakeNp,
        net_uptakeNh,

        overflowNp,
        overflowCp,
        overflowNh,
        overflowCh,

        deathp , deathh ,
        leakinessOp, leakinessIp, leakinessOh, leakinessIh, 
        ROSreleasep, ROSbreakdownh,
        respirationp, respirationh, dic_air_water_exchange,
        limICp,  exp(-omegaP*ROS),
    ]
    interm_names = [
        'Xp', 'Xh',
        'limINp',
        'limONp',
        'limICp',
        'limOCp',
        'limINh',
        'limONh',
        'limICh',
        'limOCh',

        'gross_uptakeINp',
        'gross_uptakeONp',
        'gross_uptakeICp',
        'gross_uptakeOCp',
        'gross_uptakeINh',
        'gross_uptakeONh',
        'gross_uptakeICh',
        'gross_uptakeOCh',

        'net_uptakeNp',
        'net_uptakeNh',

        'overflowNp',
        'overflowCp',
        'overflowNh',
        'overflowCh',
        'deathp' , 'deathh' ,
        'leakinessOp', 'leakinessIp', 'leakinessOh', 'leakinessIh', 
        'ROSreleasep', 'ROSbreakdownh',    
        'respirationp', 'respirationh', 'dic_air_water_exchange',
        
        'limICp',  'exp(-omegaP*ROS)',
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt,  dDONdt_ponly, dRDONdt_ponly, dDINdt_ponly, dDOCdt_ponly, dRDOCdt_ponly, dDICdt_ponly, dROSdt_ponly, dABpdt, dABhdt_ponly]
    var_list   = [ Bp,    DON,    RDON,    DIN,    DOC,    RDOC,  DIC,    ROS,    ABp,   ABh]
    var_names  = ['Bp',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'ABp', 'ABh']
    init_vars = [INIT_BP,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xp, 

        limINp,
        limONp,
        limICp,
        limOCp,

        gross_uptakeINp,
        gross_uptakeONp,
        gross_uptakeICp,
        gross_uptakeOCp,

        net_uptakeNp,

        overflowNp,
        overflowCp,

        deathp ,
        leakinessOp, leakinessIp, 
        ROSreleasep, 
        respirationp, dic_air_water_exchange,

    ]
    interm_names = [
        'Xp', 
        'limINp',
        'limONp',
        'limICp',
        'limOCp',

        'gross_uptakeINp',
        'gross_uptakeONp',
        'gross_uptakeICp',
        'gross_uptakeOCp',
        'net_uptakeNp',

        'overflowNp',
        'overflowCp',
        'deathp' , 
        'leakinessOp', 'leakinessIp', 
        'ROSreleasep',  
        'respirationp', 'dic_air_water_exchange',
        
    ]
    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals_str=param_vals):
    sfunc_list = [dBhdt, dDONdt_honly, dRDONdt_honly, dDINdt_honly, dDOCdt_honly, dRDOCdt_honly,dDICdt_honly, dROSdt_honly, dABpdt_honly, dABhdt]
    var_list   = [Bh,    DON,    RDON,    DIN,    DOC,    RDOC,   DIC,    ROS,    ABp,   ABh]
    var_names  = ['Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'ABp', 'ABh']
    init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS,INIT_SP,INIT_SH]
    param_vals = {symbols(k) : v for k,v in param_vals_str.items()}

    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xh,

        limINh,
        limONh,
        limICh,
        limOCh,

        gross_uptakeINh,
        gross_uptakeONh,
        gross_uptakeICh,
        gross_uptakeOCh,

        net_uptakeNh,

        overflowNh,
        overflowCh,

         deathh ,
        leakinessOh, leakinessIh, 
        ROSbreakdownh,
        respirationh, dic_air_water_exchange,
    ]
    interm_names = [
        'Xh',
        'limINh',
        'limONh',
        'limICh',
        'limOCh',

        'gross_uptakeINh',
        'gross_uptakeONh',
        'gross_uptakeICh',
        'gross_uptakeOCh',

        'net_uptakeNh',

        'overflowNh',
        'overflowCh',
        'deathh' ,
        'leakinessOh', 'leakinessIh', 
        'ROSbreakdownh',    
        'respirationh', 'dic_air_water_exchange',
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func


def print_dydt0(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    for i,j, k in zip(var_names, dydt0, init_vars):
        print(f'd{i}/dt = {j:.2e}, init {i} = {k:.2e}, newval = {k+j:.2e}')


def print_intermediate0(intermediate_func, interm_names, init_vars):
    for i,j in zip(interm_names, intermediate_func(*init_vars)):
        print(f'{i:<4} = {j:.2e}')


def biomass_diff0(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    V = dict(zip(var_names, dydt0))
    print (f"dBp/dt + dBh/dt + dDON/dt + dRDON/dt + dDIN/dt = { V['Bp'] + V['Bh'] + V['DON'] + V['RDON'] + V['DIN'] }")

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
    tend = days*seconds_in_day
    if t_eval is None: 
        t_eval = np.arange(tstart, tend, 3600*4)
    sol = solve_ivp(
        fun=calc_dydt, y0=init_vars,
        t_span=[tstart, tend], t_eval=t_eval, max_step=100, first_step=1)
    #print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval={t_eval})')
    #print(sol.message)
    return sol


def solver2df_ivp(sol, var_names, interm_names, intermediate_func):
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
    mdf = df.melt(id_vars=['t', 'day'])
    return df, mdf


def run_solver_ode(calc_dydt, init_vars, days, t_eval):
    tstart = 0
    tend = days*seconds_in_day
    if t_eval is None: 
        t_eval = np.arange(tstart, tend, 3600*4)
    y = odeint(
        func=calc_dydt, y0=init_vars,
        t=t_eval, hmax=100, h0=1, tfirst=True)
    #print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval={t_eval})')
    #print(sol.message)
    return t_eval, y
    #return solver2df_ode(sol, t_eval, var_names, interm_names, intermediate_func)

def solver2df_ode(sol, var_names, interm_names, intermediate_func):
    d = dict(zip(var_names, sol[1].T))
    d['t'] = sol[0]
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
    mdf = df.melt(id_vars=['t', 'day'])
    return df, mdf


def run_solver(calc_dydt, init_vars, days=140, t_eval=None):
    tstart = time.process_time()
    sol = run_solver_ode(calc_dydt, init_vars, days, t_eval=t_eval)
    tend =  time.process_time()
    print ('simulation time', tend - tstart)
    return sol
    
def solver2df(sol, var_names, interm_names, intermediate_func):
    return solver2df_ode(sol, var_names, interm_names, intermediate_func)

def get_t_eval(maxday, step = 3600*4, ref_times = None):
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


def _rmse(refdf, df, refcol, col):
    smallrefdf = refdf.dropna(subset=[refcol])
    ref_t = np.rint(smallrefdf['t'])
    tdf = df.loc[df.t.isin(ref_t)]
    return mean_squared_error(tdf[col], smallrefdf[refcol])
   

def run_with_params_json(json_fpath, days, refdf, out_dpath, out_fprefix):
    perr = -1
    herr = -1
    new_params = json2params(param_vals, json_fpath)
    var_names, init_vars, calc_dydt, interm_names, intermediate_func = get_main_data(param_vals_str=new_params)
    if refdf is not None:
        ref_t = np.rint(refdf['t'])
        t_eval = get_t_eval(days, ref_times = ref_t)
    else:
        t_eval = get_t_eval(days)
    
    sol = run_solver(calc_dydt, init_vars, t_eval=t_eval, days=days)
    sumdf = pd.DataFrame({str(k): v for k,v in new_params.items()}, index=[0])
    sumdf['run_id'] = out_fprefix
    #sumdf['status'] = sol.status
    #if sol.status != 0:
    #    sumdf['message'] = sol.message
    #if sol.success:
    df, mdf = solver2df(sol, var_names, interm_names, intermediate_func)
    df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'), compression='gzip')

    if refdf is not None:
        perr =  _rmse(refdf, df, refcol = 'cc Bp[N]', col='Bp')
        herr =  _rmse(refdf, df, refcol = 'cc Bh[N]', col='Bh')
        sumdf['h_err'] = herr
        sumdf['p_err'] = perr
    sumdf.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'), compression='gzip')
    return perr + herr
   
def generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60):
    hash_val = str(hash(tuple(params.values())))
    run_id = f'{out_fprefix}_h{hash_val}'
    json_fpath = os.path.join(json_dpath, f'{run_id}_params.json')
    params2json(params, json_fpath)
    return run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout)


def get_params(X, params_to_update, param_vals, log_params=None): 
    new_param_vals = param_vals.copy()
    if log_params is not None:
        X = [np.exp(x) if lg else x for x,lg in zip(X, log_params)]
    
    new_param_vals.update({k : v for k,v in zip(params_to_update, X)})
    return new_param_vals

def generate_json_and_run_from_X(X, params_to_update, param_vals, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60, log_params=None):
    params = get_params(X, params_to_update, param_vals, log_params)
    return generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout)



def run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout=10*60):
    try:
        result = subprocess.run(
            [sys.executable, __file__, 
             '--json', json_fpath,
             '--ref_csv', ref_csv, 
             '--run_id', run_id,
             '--outdpath', out_dpath,
            ], 
            capture_output=True, text=True, check=True, timeout=timeout,
        )
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)

        m = re.search(r'^MSE: ([\d\.]+)$', result.stdout, flags=re.MULTILINE)
        if m is not None:
            return (float(m.groups()[0]))
        
    except subprocess.CalledProcessError as err:
        print('CalledProcessError', err.returncode)
        print("stdout:", err.stdout)
        print("stderr:", err.stderr)
                
    except subprocess.TimeoutExpired as err:
        print('TimeoutExpired', err.timeout)
        print("stdout:", err.stdout)
        print("stderr:", err.stderr)
    return 1e50

def run_chunk(param_vals, param_values, params_to_update, chunk, number_of_runs, run_id, ref_csv, json_dpath, out_dpath, timeout, skip_if_found=True, log_params=None):
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
                ref_csv, json_dpath, out_dpath, out_fprefix, timeout, log_params=log_params)
            

def run_sensitivity_per_parameter(param_vals, parameter, bound, number_of_runs, run_id, ref_csv, json_dpath, out_dpath, timeout, skip_if_found=True, log_param=False):
    if log_param:
        bound = (np.log(bound[0]), np.log(bound[1]))
    for i,v in enumerate(np.linspace(bound[0], bound[1],num=number_of_runs)):
        out_fprefix = f'{run_id}_{parameter}_{i}'
        print(out_fprefix)
        generate_json_and_run_from_X(
            [v], [parameter], param_vals, 
            ref_csv, json_dpath, out_dpath, out_fprefix, timeout, log_params=[log_param])


    
if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', default='prelim bottle.csv')
    parser.add_argument('--json', help='json with param vals', default=None)
    parser.add_argument('--maxday', help='max day of simulation', type=int, default=140)

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--run_id", help="run id", required=True)

    args = parser.parse_args()
    dpath = args.outdpath
    if dpath != '':
        os.makedirs(dpath, exist_ok=True)
    refdf = pd.read_csv(args.ref_csv)
    model_name = args.run_id

    MSE_err = run_with_params_json(args.json, args.maxday, refdf, dpath, args.run_id)
    print ('\nMSE:', MSE_err)

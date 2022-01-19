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
Bp, Bh, DOC, RDOC, DIC, DON, RDON, DIN, ROS, Sp, Sh = symbols('B_p B_h DOC RDOC DIC DON RDON DIN ROS S_p S_h')


# parameters
gammaDp, gammaDh, EOp, EIp, EOh, EIh = symbols('gamma^D_p gamma^D_h E^O_p E^I_p E^O_h E^I_h')
KSp, KSh, Esp, Esh, Msp, Msh = symbols('K^S_p K^S_h E^S_p E^S_h M^S_p M^S_h')

Op, Oh = symbols('O_p O_h')
epsilon, VTmax, KTh, omega = symbols('epsilon VTmax KT_h omega')
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

# initial concentrations
INIT_DIN = 100
INIT_DON = 20
INIT_RDON = 0
INIT_RDOC = 0
INIT_DIC = 3000
INIT_DOC = INIT_DON * R_CN
INIT_BP = 1e9 * Qp
INIT_BH = 1e10 * Qh
INIT_BH_CC = 5e9 * Qh # the actual concentration in the measurements
INIT_ROS = 0
INIT_SP = 0
INIT_SH = 0

R_P = R_CN
R_H = R_CN

# DIC exchange
h = 0.3 # height in m
Kg = 0.4968 / seconds_in_day # m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
c_sat = INIT_DIC


param_vals_with_symbols = {
    # 1/d
    Mh: 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.8,         # pro death release 
    gammaDh : 0.8,         # het death release
    
    # ratio
    Rp : R_CN,
    Rh : R_CN,
    
    # 1/d
    # EOp : 0.2 / seconds_in_day,        # move to income tax # pro organic exudation 
    # EIh : 0.2 / seconds_in_day,        # move to income tax # het inorganic exudation 
    EOp : 0.1 / seconds_in_day,        # pro organic exudation 
    EIp : 0 / seconds_in_day,          # pro inorganic exudation
    EOh : 0 / seconds_in_day,          # het organic exudation 
    EIh : 0.1 / seconds_in_day,        # het inorganic exudation 
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake

    KONp : 0.17 * pro_vol**0.27, 
    KINp : 5* 0.17 * pro_vol**0.27, # x 5 based on sensitivity
    KOCp : 0.17 * pro_vol**0.27, 
    KICp : 0.17 * pro_vol**0.27, 
    KONh : 0.17 * alt_vol**0.27, 
    KINh : 0.17 * alt_vol**0.27, # / 10 based on sensitivity
    KOCh : 5*0.17 * alt_vol**0.27, # x 5 based on sensitivity
    KICh : 0.17 * alt_vol**0.27, 
    # umol N/cell/d
    # vmax = muinfp* VmaxIp * Qp
    # 1/day  * umol/cell  * umol/cell/d # TODO - figure out units
    VmaxONp : 0.7 * 1.9e-9 / 10000 / Qp / seconds_in_day, 
    VmaxINp : 0.5 * 1.9e-9 / Qp / seconds_in_day, 
    VmaxOCp : 0.7 * 1.9e-9 / 10000 * R_P / Qp / seconds_in_day, 
    VmaxICp : 0.7 * 1.9e-9 * R_P / Qp / seconds_in_day, 
    VmaxONh : 2 * 1.9e-9 / Qh / seconds_in_day, 
    VmaxINh : 4 * 1.9e-9 / Qh / seconds_in_day, # x 3 based on sensitivity
    VmaxOCh : 2 * 1.9e-9 * R_H / Qh / seconds_in_day, 
    VmaxICh : 2 * 1.9e-9 / 10000 * R_H / Qh / seconds_in_day, 
    
    #Vmaxp : 0.8 * 1.9e-9 * pro_vol**0.67 /Qp / seconds_in_day, 
    #Vmaxh : 3 * 1.9e-9 * alt_vol**0.67 / Qh / seconds_in_day, 

    # overflow rate
    Op : 0.6, # changed based on sensitivity
    Oh : 0.6, # changed based on sensitivity
    
    # umol/cell/d
    epsilon : 1e-10 / Qp / seconds_in_day, 
    VTmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 0.01, #0.1,
    
    # Signals
    # umolN/L
    KSp : 0.17 * pro_vol**0.27 * 100, 
    KSh : 0.17 * pro_vol**0.27 * 100, 
    # umol/cell/d
    Esp : 1e-20 / Qp / seconds_in_day, 
    Esh : 1e-20 / Qh / seconds_in_day, 
    # 1/d (between -1 to 1)  / seconds_in_day
    Msp : -0.01 / seconds_in_day, 
    Msh : -0.01 / seconds_in_day,
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
    Rp : R_CN,
    Rh : R_CN,
    
    # 1/d
    # EOp : 0.2 / seconds_in_day,        # move to income tax # pro organic exudation 
    # EIh : 0.2 / seconds_in_day,        # move to income tax # het inorganic exudation 
    EOp : 0.1 / seconds_in_day,        # pro organic exudation 
    EIp : 0 / seconds_in_day,          # pro inorganic exudation
    EOh : 0 / seconds_in_day,          # het organic exudation 
    EIh : 0.1 / seconds_in_day,        # het inorganic exudation 
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

    # overflow rate
    Op : 0.6, # changed based on sensitivity
    Oh : 0.6, # changed based on sensitivity
    
    # umol/cell/d
    epsilon : 1e-10 / Qp / seconds_in_day, 
    VTmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 0.01, #0.1,
    
    # Signals
    # umolN/L
    KSp : 0.17 * pro_vol**0.27 * 100, 
    KSh : 0.17 * pro_vol**0.27 * 100, 
    # umol/cell/d
    Esp : 1e-20 / Qp / seconds_in_day, 
    Esh : 1e-20 / Qh / seconds_in_day, 
    # 1/d (between -1 to 1)  / seconds_in_day
    Msp : -0.01 / seconds_in_day, 
    Msh : -0.01 / seconds_in_day,
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
        new_param_vals[str(omega)] = 0
    if DISABLE_MECHANISMS.P_SIGNAL in mechanisms:
        new_param_vals[str(Esp)] = 0
    if DISABLE_MECHANISMS.H_SIGNAL in mechanisms:
        new_param_vals[str(Esh)] = 0

    return new_param_vals



# different model configurations
DISABLE_ROS_SIGNAL = {
    # ROS
    str(epsilon), 
    str(VTmax), 
    str(KTh), 
    str(omega), #0.1,
    
    # Signals
    # umolN/L
    str(KSp), 
    str(KSh), 
    # umol/cell/d
    str(Esp), 
    str(Esh), 
    # 1/d (between -1 to 1)  / seconds_in_day
    str(Msp), 
    str(Msh),
}

DISABLE_MIXOTROPHY = {
    str(VmaxONp), 
    str(VmaxOCp), 
}

DISABLE_EXUDATION_OVERFLOW = {
    str(Op), 
    str(Oh), 
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

# cells/L
Xp = Bp / Qp
Xh = Bh / Qh


# ratio
limINp = (DIN / (DIN + KINp))
limONp = (DON / (DON + KONp))
limICp = (DIC / (DIC + KICp))
limOCp = (DOC / (DOC + KOCp))
limINh = (DIN / (DIN + KINh))
limONh = (DON / (DON + KONh))
limICh = (DIC / (DIC + KICh))
limOCh = (DOC / (DOC + KOCh))

limSh = (Sh / (Sh + KSh))
limSp = (Sp / (Sp + KSp))


# vmax = muinfp* VmaxIp / Qp
# umol N /L 
gross_uptakeINp = VmaxINp * limINp * exp(-omega*ROS) * Bp
gross_uptakeONp = VmaxONp * limONp * exp(-omega*ROS) * Bp
gross_uptakeINh = VmaxINh * limINh * Bh
gross_uptakeONh = VmaxONh * limONh * Bh
# umol C /L
gross_uptakeICp = VmaxICp * limICp * exp(-omega*ROS) * Bp 
gross_uptakeOCp = VmaxOCp * limOCp * exp(-omega*ROS) * Bp 
gross_uptakeICh = VmaxICh * limICh * Bh 
gross_uptakeOCh = VmaxOCh * limOCh * Bh 

# umol N / L
actual_uptakeNp = Min(gross_uptakeINp + gross_uptakeONp, 
                     (gross_uptakeICp + gross_uptakeOCp) / Rp)

actual_uptakeNh = Min(gross_uptakeINh + gross_uptakeONh, 
                     (gross_uptakeICh + gross_uptakeOCh) / Rh)

# ratio
IOuptakeRateNp = gross_uptakeINp / (gross_uptakeINp + gross_uptakeONp)
IOuptakeRateCp = gross_uptakeICp / (gross_uptakeICp + gross_uptakeOCp)
IOuptakeRateNh = gross_uptakeINh / (gross_uptakeINh + gross_uptakeONh)
IOuptakeRateCh = gross_uptakeICh / (gross_uptakeICh + gross_uptakeOCh)

# umol N / L
overflowNp = gross_uptakeINp + gross_uptakeONp - actual_uptakeNp 
overflowNh = gross_uptakeINh + gross_uptakeONh - actual_uptakeNh
# umol C /L
overflowCp = gross_uptakeICp + gross_uptakeOCp - actual_uptakeNp * Rp
overflowCh = gross_uptakeICh + gross_uptakeOCh - actual_uptakeNh * Rh

# umol N / L
overflowINp =                   (1 - Op) * overflowNp * IOuptakeRateNp
overflowONp = Op * overflowNp + (1 - Op) * overflowNp * (1 - IOuptakeRateNp)

overflowICp =                   (1 - Op) * overflowCp * IOuptakeRateCp
overflowOCp = Op * overflowCp + (1 - Op) * overflowCp * (1 - IOuptakeRateCp)

overflowINh = Oh * overflowNh + (1 - Oh) * overflowNh * IOuptakeRateNh
overflowONh =                   (1 - Oh) * overflowNh * (1 - IOuptakeRateNh)
overflowICh = Oh * overflowCh + (1 - Oh) * overflowCh * IOuptakeRateCh
overflowOCh =                   (1 - Oh) * overflowCh * (1 - IOuptakeRateCh)

respirationp =  bp* actual_uptakeNp + Bp * r0p
respirationh =  bh* actual_uptakeNh + Bh * r0h
dic_uptake   =  - (DIC - c_sat) / tau




# M = M / Q

# death = M * X = M * B /Q = M / Q * B

death_ratep = Min(Max(Mp - Msp*limSh, 0), 1 / seconds_in_day)
death_rateh = Min(Max(Mh - Msh*limSp, 0), 1 / seconds_in_day)

deathp = death_ratep * Bp #* Xp
deathh = death_rateh * Bh #* Xh

exudationOp = EOp * Bp
exudationIp = EIp * Bp
exudationOh = EOh * Bh
exudationIh = EIh * Bh

# epsilon = epsilon / Q
# VTMax = VTmax / Q
Treleasep = epsilon * Bp
Tbreakdownh = VTmax * ROS / (ROS + KTh) * Bh

# signal
Sreleasep = Esp * Bp
Sreleaseh = Esh * Bh


# final equation - coculture
dBpdt = actual_uptakeNp - deathp - exudationOp - exudationIp - respirationp - Sreleasep
dBhdt = actual_uptakeNh - deathh - exudationOh - exudationIh - respirationh - Sreleaseh
dDONdt = deathp * gammaDp + deathh * gammaDh + exudationOp +  exudationOh - gross_uptakeONp -  gross_uptakeONh + overflowONp + overflowONh 
dDOCdt = deathp * gammaDp * Rp + deathh * gammaDh * Rh + exudationOp *Rp +  exudationOh * Rh - gross_uptakeOCp -  gross_uptakeOCh + overflowOCp + overflowOCh
dRDONdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dRDOCdt = deathp * (1 - gammaDp) * Rp + deathh * (1 - gammaDh) * Rh

dDINdt = exudationIp +  exudationIh - gross_uptakeINp -  gross_uptakeINh  + overflowINp + overflowINh + respirationh + respirationp 
dDICdt = exudationIp *Rp +  exudationIh* Rh - gross_uptakeICp -  gross_uptakeICh  + overflowICp + overflowICh + respirationh* Rh + respirationp * Rp + dic_uptake

dROSdt = Treleasep - Tbreakdownh 
dSpdt = Sreleasep
dShdt = Sreleaseh



# PRO only model
dDONdt_ponly = deathp * gammaDp + exudationOp - gross_uptakeONp  + overflowONp  
dDOCdt_ponly = deathp * gammaDp * Rp + exudationOp *Rp - gross_uptakeOCp + overflowOCp 
dRDONdt_ponly = deathp * (1 - gammaDp) 
dRDOCdt_ponly = deathp * (1 - gammaDp) * Rp
dDINdt_ponly = exudationIp - gross_uptakeINp + overflowINp+ respirationp
dDICdt_ponly = exudationIp *Rp - gross_uptakeICp + overflowICp + respirationp* Rp + dic_uptake
dROSdt_ponly = Treleasep  
dShdt_ponly = Integer(0)

# HET only model
dDONdt_honly = deathh * gammaDh +  exudationOh -  gross_uptakeONh + overflowONh
dDOCdt_honly = deathh * gammaDh * Rh +  exudationOh * Rh - gross_uptakeOCh + overflowOCh
dRDONdt_honly = deathh * (1 - gammaDh)
dRDOCdt_honly = deathh * (1 - gammaDh) * Rh
dDINdt_honly = exudationIh -  gross_uptakeINh  + overflowINh + respirationh
dDICdt_honly = exudationIh* Rh -  gross_uptakeICh  + overflowICh + respirationh* Rh + dic_uptake
dROSdt_honly = - Tbreakdownh 
dSpdt_honly = Integer(0)


def print_equations():
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'Sp', 'Sh']
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dSpdt, dShdt]
    for n,f in zip(var_names, sfunc_list):
        print(f'd{n}/dt')
        display(f)

def get_main_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt, dSpdt, dShdt]
    var_list   = [ Bp,    Bh,    DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS,    Sp,    Sh]
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS',   'Sp',  'Sh']
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

        actual_uptakeNp,
        actual_uptakeNh,

        overflowNp,
        overflowCp,
        overflowNh,
        overflowCh,

        deathp , deathh ,
        exudationOp, exudationIp, exudationOh, exudationIh, 
        Treleasep, Tbreakdownh,
        respirationp, respirationh, dic_uptake,
        
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

        'actual_uptakeNp',
        'actual_uptakeNh',

        'overflowNp',
        'overflowCp',
        'overflowNh',
        'overflowCh',
        'deathp' , 'deathh' ,
        'exudationOp', 'exudationIp', 'exudationOh', 'exudationIh', 
        'Treleasep', 'Tbreakdownh',    
        'respirationp', 'respirationh', 'dic_uptake',
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt,  dDONdt_ponly, dRDONdt_ponly, dDINdt_ponly, dDOCdt_ponly, dRDOCdt_ponly, dDICdt_ponly, dROSdt_ponly, dSpdt, dShdt_ponly]
    var_list   = [ Bp,    DON,    RDON,    DIN,    DOC,    RDOC,  DIC,    ROS,    Sp,   Sh]
    var_names  = ['Bp',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'Sp', 'Sh']
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

        actual_uptakeNp,

        overflowNp,
        overflowCp,

        deathp ,
        exudationOp, exudationIp, 
        Treleasep, 
        respirationp, dic_uptake,

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
        'actual_uptakeNp',

        'overflowNp',
        'overflowCp',
        'deathp' , 
        'exudationOp', 'exudationIp', 
        'Treleasep',  
        'respirationp', 'dic_uptake',
        
    ]
    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals_str=param_vals):
    sfunc_list = [dBhdt, dDONdt_honly, dRDONdt_honly, dDINdt_honly, dDOCdt_honly, dRDOCdt_honly,dDICdt_honly, dROSdt_honly, dSpdt_honly, dShdt]
    var_list   = [Bh,    DON,    RDON,    DIN,    DOC,    RDOC,   DIC,    ROS,    Sp,   Sh]
    var_names  = ['Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS', 'Sp', 'Sh']
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

        actual_uptakeNh,

        overflowNh,
        overflowCh,

         deathh ,
        exudationOh, exudationIh, 
        Tbreakdownh,
        respirationh, dic_uptake,
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

        'actual_uptakeNh',

        'overflowNh',
        'overflowCh',
        'deathh' ,
        'exudationOh', 'exudationIh', 
        'Tbreakdownh',    
        'respirationh', 'dic_uptake',
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
    if 'Sp' in df.columns:
        df['Sp[C]'] = df['Sp']*param_vals[str(Rp)]
    if 'Sh' in df.columns:
        df['Sh[C]'] = df['Sh']*param_vals[str(Rh)]
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
    if 'Sp' in df.columns:
        df['Sp[C]'] = df['Sp']*param_vals[str(Rp)]
    if 'Sh' in df.columns:
        df['Sh[C]'] = df['Sh']*param_vals[str(Rh)]
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

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
init_printing(use_unicode=True)
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error
import json
import subprocess
import sys
import re

# variables
Bp, Bh, DOC, RDOC, DIC, DON, RDON, DIN, ROS = symbols('B_p B_h DOC RDOC DIC DON RDON DIN ROS')

# parameters
gammaDp, gammaDh, EOp, EIp, EOh, EIh = symbols('gamma^D_p gamma^D_h E^O_p E^I_p E^O_h E^I_h')

Op, Oh = symbols('O_p O_h')
epsilon, VTmax, KTh, omega = symbols('epsilon VTmax KT_h omega')
Mp, Mh = symbols('M_p M_h')
Rp, Rh = symbols('R_p R_h')

KONp, KINp, KOCp, KICp, KONh, KINh, KOCh, KICh = symbols('K^ON_p K^IN_p K^OC_p K^IC_p K^ON_h K^IN_h K^OC_h K^IC_h')
VmaxONp, VmaxINp, VmaxOCp, VmaxICp, VmaxONh, VmaxINh, VmaxOCh, VmaxICh = symbols('Vmax^ON_p Vmax^IN_p Vmax^OC_p Vmax^IC_p Vmax^ON_h Vmax^IN_h Vmax^OC_h Vmax^IC_h')


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

R_P = R_CN
R_H = R_CN

param_vals_with_symbols = {
    # 1/d
    Mh: 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.5,         # pro death release 
    gammaDh : 0.5,         # het death release
    
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

    KONp : 0.17 * pro_vol**0.27 * 10, 
    KINp : 0.17 * pro_vol**0.27 * 10, 
    KOCp : 0.17 * pro_vol**0.27 * 10, 
    KICp : 0.17 * pro_vol**0.27 * 10, 
    KONh : 0.17 * alt_vol**0.27 * 10, 
    KINh : 0.17 * alt_vol**0.27 * 10, 
    KOCh : 0.17 * alt_vol**0.27 * 10, 
    KICh : 0.17 * alt_vol**0.27 * 10, 
    # umol N/cell/d
    # vmax = muinfp* VmaxIp * Qp
    # 1/day  * umol/cell  * umol/cell/d # TODO - figure out units
    VmaxONp : 0.7 * 1.9e-9 / 10000 / Qp / seconds_in_day, 
    VmaxINp : 0.5 * 1.9e-9 / Qp / seconds_in_day, 
    VmaxOCp : 0.7 * 1.9e-9 / 10000 * R_P / Qp / seconds_in_day, 
    VmaxICp : 0.7 * 1.9e-9 * R_P / Qp / seconds_in_day, 
    VmaxONh : 2 * 1.9e-9 / Qh / seconds_in_day, 
    VmaxINh : 2 * 1.9e-9 / Qh / seconds_in_day, 
    VmaxOCh : 2 * 1.9e-9 * R_H / Qh / seconds_in_day, 
    VmaxICh : 2 * 1.9e-9 / 10000 * R_H / Qh / seconds_in_day, 
    
    #Vmaxp : 0.8 * 1.9e-9 * pro_vol**0.67 /Qp / seconds_in_day, 
    #Vmaxh : 3 * 1.9e-9 * alt_vol**0.67 / Qh / seconds_in_day, 

    # overflow rate
    Op : 1,
    Oh : 1,
    
    # umol/cell/d
    epsilon : 1e-10 / Qp / seconds_in_day, 
    VTmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 0.1,
    
}
param_vals = {str(k) : v for k,v in param_vals_with_symbols.items()}


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
    #PROALLELO = auto()
    #HETALLELO = auto()


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
    #if PRO_ALLELO in mechanisms:
        # TODO
        # new_param_vals[VmaxINh] = new_param_vals[VmaxINh] * 1e-2
    #    pass
    #if HET_ALLELO in mechanisms:
        # TODO
        # new_param_vals[VmaxINh] = new_param_vals[VmaxINh] * 1e-2
    #    pass

    return new_param_vals



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




# M = M / Q

# death = M * X = M * B /Q = M / Q * B

deathp = Mp * Bp #* Xp
deathh = Mh * Bh #* Xh
exudationOp = EOp * Bp
exudationIp = EIp * Bp
exudationOh = EOh * Bh
exudationIh = EIh * Bh

# epsilon = epsilon / Q
# VTMax = VTmax / Q
Treleasep = epsilon * Bp
Tbreakdownh = VTmax * ROS / (ROS + KTh) * Bh

dBpdt = actual_uptakeNp - deathp - exudationOp - exudationIp
dBhdt = actual_uptakeNh - deathh - exudationOh - exudationIh
dDONdt = deathp * gammaDp + deathh * gammaDh + exudationOp +  exudationOh - gross_uptakeONp -  gross_uptakeONh + overflowONp + overflowONh 
dDOCdt = deathp * gammaDp * Rp + deathh * gammaDh * Rh + exudationOp *Rp +  exudationOh * Rh - gross_uptakeOCp -  gross_uptakeOCh + overflowOCp + overflowOCh
dRDONdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dRDOCdt = deathp * (1 - gammaDp) * Rp + deathh * (1 - gammaDh) * Rh

dDINdt = exudationIp +  exudationIh - gross_uptakeINp -  gross_uptakeINh  + overflowINp + overflowINh 
dDICdt = exudationIp *Rp +  exudationIh* Rh - gross_uptakeICp -  gross_uptakeICh  + overflowICp + overflowICh

dROSdt = Treleasep - Tbreakdownh 

# PRO only model
dDONdt_ponly = deathp * gammaDp + exudationOp - gross_uptakeONp  + overflowONp  
dDOCdt_ponly = deathp * gammaDp * Rp + exudationOp *Rp - gross_uptakeOCp + overflowOCp 
dRDONdt_ponly = deathp * (1 - gammaDp) 
dRDOCdt_ponly = deathp * (1 - gammaDp) * Rp
# todo if needed: dRDOCdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDINdt_ponly = exudationIp - gross_uptakeINp + overflowINp
dDICdt_ponly = exudationIp *Rp - gross_uptakeICp + overflowICp
dROSdt_ponly = Treleasep  

# HET only model
dDONdt_honly = deathh * gammaDh +  exudationOh -  gross_uptakeONh + overflowONh
dDOCdt_honly = deathh * gammaDh * Rh +  exudationOh * Rh - gross_uptakeOCh + overflowOCh
dRDONdt_honly = deathh * (1 - gammaDh)
dRDOCdt_honly = deathh * (1 - gammaDh) * Rh
# todo if needed: dRDOCdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDINdt_honly = exudationIh -  gross_uptakeINh  + overflowINh 
dDICdt_honly = exudationIh* Rh -  gross_uptakeICh  + overflowICh
dROSdt_honly = - Tbreakdownh 


def print_equations():
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt]
    for n,f in zip(var_names, sfunc_list):
        print(f'd{n}/dt')
        display(f)

def get_main_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dRDOCdt, dDICdt, dROSdt]
    var_list   = [ Bp,    Bh,    DON,    RDON,    DIN,    DOC,   RDOC,    DIC,    ROS]
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    init_vars = [INIT_BP,INIT_BH_CC,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
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
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_data(param_vals_str=param_vals):
    sfunc_list = [dBpdt,  dDONdt_ponly, dRDONdt_ponly, dDINdt_ponly, dDOCdt_ponly, dRDOCdt_ponly, dDICdt_ponly, dROSdt_ponly]
    var_list   = [ Bp,    DON,    RDON,    DIN,    DOC,    RDOC,  DIC,    ROS]
    var_names  = ['Bp',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    init_vars = [INIT_BP,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC,INIT_DIC,INIT_ROS]
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
    ]
    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals_str=param_vals):
    sfunc_list = [dBhdt, dDONdt_honly, dRDONdt_honly, dDINdt_honly, dDOCdt_honly, dRDOCdt_honly,dDICdt_honly, dROSdt_honly]
    var_list   = [Bh,    DON,    RDON,    DIN,    DOC,    RDOC,   DIC,    ROS]
    var_names  = ['Bh',  'DON',  'RDON',  'DIN',  'DOC',  'RDOC', 'DIC',  'ROS']
    init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_RDOC, INIT_DIC,INIT_ROS]
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



def run_solver(calc_dydt, init_vars, days=63, t_eval=None):
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


def solver2df(sol, var_names, interm_names, intermediate_func):
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
    mdf = df.melt(id_vars=['t', 'day'])
    return df, mdf


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




def run_with_params_json(json_fpath, days, refdf, out_dpath, out_fprefix):
    perr = -1
    herr = -1
    new_params = json2params(param_vals, json_fpath)
    var_names, init_vars, calc_dydt, interm_names, intermediate_func = get_main_data(param_vals_str=new_params)
    t_eval=None
    if refdf is not None:
        t_eval=refdf.t.values
    sol = run_solver(calc_dydt, init_vars, t_eval=t_eval, days=days)
    sumdf = pd.DataFrame({str(k): v for k,v in new_params.items()}, index=[0])
    sumdf['run_id'] = out_fprefix
    sumdf['status'] = sol.status
    if sol.status != 0:
        sumdf['message'] = sol.message
    if sol.success:
        df, mdf = solver2df(sol, var_names, interm_names, intermediate_func)
        df.to_csv(os.path.join(out_dpath, f'{out_fprefix}_df.csv.gz'))

        if refdf is not None:
            perr = mean_squared_error(df.Bp, refdf['cc Bp[N]'])
            herr = mean_squared_error(df.Bh, refdf['cc Bh[N]'])
            sumdf['h_err'] = herr
            sumdf['p_err'] = perr
    sumdf.to_csv(os.path.join(out_dpath, f'{out_fprefix}_sum.csv.gz'))
    return perr + herr
   
def generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60):
    hash_val = str(hash(tuple(params.values())))
    run_id = f'{out_fprefix}_h{hash_val}'
    json_fpath = os.path.join(json_dpath, f'{run_id}_params.json')
    params2json(params, json_fpath)
    return run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout)

def get_params(X, params_to_update, param_vals): 
    new_param_vals = param_vals.copy()
    new_param_vals.update({k : v for k,v in zip(params_to_update, X)})
    return new_param_vals

def generate_json_and_run_from_X(X, params_to_update, param_vals, ref_csv, json_dpath, out_dpath, out_fprefix, timeout=10*60):
    params = get_params(X, params_to_update, param_vals)
    return generate_json_and_run(params, ref_csv, json_dpath, out_dpath, out_fprefix, timeout)



def run_with_timout(json_fpath, ref_csv, out_dpath, run_id, timeout=10*60):
    print('file', __file__)
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


if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--ref_csv', help='reference CSV', default='prelim bottle.csv')
    parser.add_argument('--json', help='json with param vals', default=None)
    parser.add_argument('--maxday', help='max day of simulation', type=int, default=63)

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

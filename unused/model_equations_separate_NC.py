#!/usr/bin/env python
# coding: utf-8

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

# variables
Bp, Bh, DOC, DIC, DON, RDON, DIN, T = symbols('B_p B_h DOC DIC DON RDON DIN T')

# parameters
gammaDp, gammaDh, EOp, EIp, EOh, EIh = symbols('gamma^D_p gamma^D_h E^O_p E^I_p E^O_h E^I_h')
Vmaxp, Vmaxh = symbols('Vmax_p Vmax_h')
muinfp, muinfh, epsilon, VTmax, KTh, omega = symbols('mu_inf_p mu_inf_h epsilon VTmax KT_h omega')
Mp, Mh = symbols('M_p M_h')
Rp, Rh = symbols('R_p R_h')

KONp, KINp, KOCp, KICp, KONh, KINh, KOCh, KICh = symbols('K^ON_p K^IN_p K^OC_p K^IC_p K^ON_h K^IN_h K^OC_h K^IC_h')


# Redfield ratio
R_CN = 6.625

# parameter values
pro_radius = 0.3628;  # "MED4" = 9312
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
INIT_DIC = 3000
INIT_DOC = INIT_DON * R_CN
INIT_BP = 1e9 * Qp
INIT_BH = 1e10 * Qh
INIT_T = 0

param_vals = {
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
    EOp : 0.2 / seconds_in_day,        # pro organic exudation 
    EIp : 0 / seconds_in_day,          # pro inorganic exudation
    EOh : 0 / seconds_in_day,          # het organic exudation 
    EIh : 0.2 / seconds_in_day,        # het inorganic exudation 
    # TODO - change k for organic/inorganic
    # umol/l
    # K = 0 --> 1
    # K = N --> ½
    # K >> N --> >> 0 --> no uptake
    # K << N --> >> 1 --> max uptake

    KONp : 1000, # no uptake # 0.17 * pro_vol**0.27, 
    KINp : 0.17 * pro_vol**0.27, 
    KOCp : 10000, # no uptake # 0.17 * pro_vol**0.27, 
    KICp : 0.17 * pro_vol**0.27, 
    KONh : 0.17 * alt_vol**0.27, 
    KINh : 0.17 * alt_vol**0.27, 
    KOCh : 0.17 * alt_vol**0.27, 
    KICh : 10000, # no uptake # 0.17 * alt_vol**0.27, 
    # umol N/cell/d
    # vmax = muinfp* VmaxIp * Qp
    # 1/day  * umol/cell  * umol/cell/d # TODO - figure out units
    Vmaxp : 0.6 * 1.9e-9 / Qp / seconds_in_day, 
    Vmaxh : 5 * 1.9e-9 / Qh / seconds_in_day, 
    
    #Vmaxp : 0.8 * 1.9e-9 * pro_vol**0.67 /Qp / seconds_in_day, 
    #Vmaxh : 3 * 1.9e-9 * alt_vol**0.67 / Qh / seconds_in_day, 
    
    # umol/cell/d
    epsilon : 1e-10 / Qp / seconds_in_day, 
    VTmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 0,
}

def print_params(param_vals=param_vals):
    for i in param_vals:
        print(i, f' = {param_vals[i]:.2e}')

# functions

Xp = Bp / Qp
Xh = Bh / Qh


# vmax = muinfp* VmaxIp / Qp
limINp = (DIN / (DIN + KINp))
limONp = (DON / (DON + KONp))
limICp = (DIC / (DIC + KICp))
limOCp = (DOC / (DOC + KOCp))
limINh = (DIN / (DIN + KINh))
limONh = (DON / (DON + KONh))
limICh = (DIC / (DIC + KICh))
limOCh = (DOC / (DOC + KOCh))

min_limp = Min(limINp + limONp, limICp + limOCp)
min_limh = Min(limINh + limONh, limICh + limOCh)

gross_uptakeINp = Vmaxp * limINp * exp(-omega*T) * Bp
gross_uptakeONp = Vmaxp * limONp * exp(-omega*T) * Bp
gross_uptakeICp = Vmaxp * limICp * exp(-omega*T) * Bp * Rp
gross_uptakeOCp = Vmaxp * limOCp * exp(-omega*T) * Bp * Rp
gross_uptakeINh = Vmaxh * limINh * Bh
gross_uptakeONh = Vmaxh * limONh * Bh
gross_uptakeICh = Vmaxh * limICh * Bh * Rh
gross_uptakeOCh = Vmaxh * limOCh * Bh * Rh

actual_uptakep = Vmaxp * min_limp * exp(-omega*T) * Bp
actual_uptakeh = Vmaxh * min_limh * Bh

overflowNp = gross_uptakeINp + gross_uptakeONp - actual_uptakep 
overflowCp = gross_uptakeICp + gross_uptakeOCp - actual_uptakep * Rp
overflowNh = gross_uptakeINh + gross_uptakeONh - actual_uptakeh 
overflowCh = gross_uptakeICh + gross_uptakeOCh - actual_uptakeh * Rh

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
Tbreakdownh = VTmax * T / (T + KTh) * Bh

dBpdt = actual_uptakep - deathp - exudationOp - exudationIp
dBhdt = actual_uptakeh - deathh - exudationOh - exudationIh
dDONdt = deathp * gammaDp + deathh * gammaDh + exudationOp +  exudationOh - gross_uptakeONp -  gross_uptakeONh + overflowNp  
dDOCdt = deathp * gammaDp * Rp + deathh * gammaDh * Rh + exudationOp *Rp +  exudationOh * Rh - gross_uptakeOCp -  gross_uptakeOCh + overflowCp 
dRDONdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
# todo if needed: dRDOCdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDINdt = exudationIp +  exudationIh - gross_uptakeINp -  gross_uptakeINh  + overflowNh 
dDICdt = exudationIp *Rp +  exudationIh* Rh - gross_uptakeICp -  gross_uptakeICh  + overflowCh

dTdt = Treleasep - Tbreakdownh 

dDONdt_ponly = deathp * gammaDp + exudationOp - gross_uptakeONp  + overflowNp  
dDOCdt_ponly = deathp * gammaDp * Rp + exudationOp *Rp - gross_uptakeOCp + overflowCp 
dRDONdt_ponly = deathp * (1 - gammaDp) 
# todo if needed: dRDOCdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDINdt_ponly = exudationIp - gross_uptakeINp 
dDICdt_ponly = exudationIp *Rp - gross_uptakeICp 
dTdt_ponly = Treleasep  

dDONdt_honly = deathh * gammaDh +  exudationOh -  gross_uptakeONh 
dDOCdt_honly = deathh * gammaDh * Rh +  exudationOh * Rh - gross_uptakeOCh 
dRDONdt_honly = deathh * (1 - gammaDh)
# todo if needed: dRDOCdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDINdt_honly = exudationIh -  gross_uptakeINh  + overflowNh 
dDICdt_honly = exudationIh* Rh -  gross_uptakeICh  + overflowCh
dTdt_honly = - Tbreakdownh 


def print_equations():
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'DIC',  'T']
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dDICdt, dTdt]
    for n,f in zip(var_names, sfunc_list):
        print(f'd{n}/dt')
        display(f)

def get_main_data(param_vals=param_vals):
    sfunc_list = [dBpdt, dBhdt, dDONdt, dRDONdt, dDINdt, dDOCdt, dDICdt, dTdt]
    var_list   = [ Bp,    Bh,    DON,    RDON,    DIN,    DOC,    DIC,    T]
    var_names  = ['Bp',  'Bh',  'DON',  'RDON',  'DIN',  'DOC',  'DIC',  'T']
    init_vars = [INIT_BP,INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_DIC,INIT_T]
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

        min_limp,
        min_limh,

        gross_uptakeINp,
        gross_uptakeONp,
        gross_uptakeICp,
        gross_uptakeOCp,
        gross_uptakeINh,
        gross_uptakeONh,
        gross_uptakeICh,
        gross_uptakeOCh,

        actual_uptakep,
        actual_uptakeh,

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

        'min_limp',
        'min_limh',

        'gross_uptakeINp',
        'gross_uptakeONp',
        'gross_uptakeICp',
        'gross_uptakeOCp',
        'gross_uptakeINh',
        'gross_uptakeONh',
        'gross_uptakeICh',
        'gross_uptakeOCh',

        'actual_uptakep',
        'actual_uptakeh',

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

def get_ponly_data(param_vals=param_vals):
    sfunc_list = [dBpdt,  dDONdt_ponly, dRDONdt_ponly, dDINdt_ponly, dDOCdt_ponly, dDICdt_ponly, dTdt_ponly]
    var_list   = [ Bp,    DON,    RDON,    DIN,    DOC,    DIC,    T]
    var_names  = ['Bp',  'DON',  'RDON',  'DIN',  'DOC',  'DIC',  'T']
    init_vars = [INIT_BP,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_DIC,INIT_T]
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xp, 

        limINp,
        limONp,
        limICp,
        limOCp,

        min_limp,

        gross_uptakeINp,
        gross_uptakeONp,
        gross_uptakeICp,
        gross_uptakeOCp,

        actual_uptakep,

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

        'min_limp',

        'gross_uptakeINp',
        'gross_uptakeONp',
        'gross_uptakeICp',
        'gross_uptakeOCp',
        'actual_uptakep',

        'overflowNp',
        'overflowCp',
        'deathp' , 
        'exudationOp', 'exudationIp', 
        'Treleasep',  
    ]
    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify(var_list, interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals=param_vals):
    sfunc_list = [dBhdt, dDONdt_honly, dRDONdt_honly, dDINdt_honly, dDOCdt_honly, dDICdt_honly, dTdt_honly]
    var_list   = [Bh,    DON,    RDON,    DIN,    DOC,    DIC,    T]
    var_names  = ['Bh',  'DON',  'RDON',  'DIN',  'DOC',  'DIC',  'T']
    init_vars = [INIT_BH,INIT_DON,INIT_RDON,INIT_DIN,INIT_DOC,INIT_DIC,INIT_T]
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify(var_list, subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xh,

        limINh,
        limONh,
        limICh,
        limOCh,

        min_limh,

        gross_uptakeINh,
        gross_uptakeONh,
        gross_uptakeICh,
        gross_uptakeOCh,

        actual_uptakeh,

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

        'min_limh',

        'gross_uptakeINh',
        'gross_uptakeONh',
        'gross_uptakeICh',
        'gross_uptakeOCh',

        'actual_uptakeh',

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



def run_solver(calc_dydt, init_vars, days=62, t_eval=None):
    tstart = 0
    tend = days*seconds_in_day
    if t_eval is None: 
        t_eval = np.arange(tstart, tend, 3600*4)
    sol = solve_ivp(
        fun=calc_dydt, y0=init_vars,
        t_span=[tstart, tend], t_eval=t_eval, max_step=100, first_step=1)
    print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval={t_eval})')
    print(sol.message)
    return sol


def solver2df(sol, var_names, interm_names, intermediate_func):
    d = dict(zip(var_names, sol.y))
    d['t'] = sol.t
    df = pd.DataFrame(data=d)
    df['day'] = df['t']/seconds_in_day
    df[interm_names] = df[var_names].apply(lambda x : intermediate_func(*x), axis=1, 
                                           result_type='expand')
    mdf = df.melt(id_vars=['t', 'day'])
    return df, mdf


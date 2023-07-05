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
Bp, Bh, DOM, RDOM, DIM, T = symbols('B_p B_h DOM RDOM DIM T')

# parameters
gammaDp, gammaDh, EOp, EIp, EOh, EIh = symbols('gamma^D_p gamma^D_h E^O_p E^I_p E^O_h E^I_h')
Qmaxp, Qminp, Qmaxh, Qminh, KOp, KIp, KOh, KIh, VmaxOp, VmaxIp, VmaxOh, VmaxIh = symbols('Qmax_p Qmin_p Qmax_h Qmin_h KO_p KI_p KO_h KI_h VmaxO_p VmaxI_p VmaxO_h VmaxI_h')
muinfp, muinfh, epsilon, VTmax, KTh, omega = symbols('mu_inf_p mu_inf_h epsilon VTmax KT_h omega')
Mp, Mh = symbols('M_p M_h')

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


param_vals = {
    # 1/d
    Mh: 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.8,         # pro death release 
    gammaDh : 0.8,         # het death release
    # 1/d
    EOp : 0.2 / seconds_in_day,        # pro organic exudation 
    EIp : 0 / seconds_in_day,          # pro inorganic exudation
    EOh : 0 / seconds_in_day,          # het organic exudation 
    EIh : 0.2 / seconds_in_day,        # het inorganic exudation 
    # TODO - change k for organic/inorganic
    # umol/l
    KOp : 0.17 * pro_vol**0.27, 
    KIp : 0.17 * pro_vol**0.27, 
    KOh : 0.17 * alt_vol**0.27, 
    KIh : 0.17 * alt_vol**0.27, 
    # umol/cell/d
    # vmax = muinfp* VmaxIp / Qp
    VmaxOp : 0, 
    VmaxIp : 0.6 * 1.9e-9 / Qp / seconds_in_day, 
    VmaxOh : 1.5 * 1.9e-9* pro_alt_vol_ratio**0.67 / Qh / seconds_in_day, 
    VmaxIh : 1.5 * 1.9e-9 / 10* pro_alt_vol_ratio**0.67 / Qh / seconds_in_day,
    # umol/cell/d
    epsilon : 1e-10 / Qp / seconds_in_day, 
    VTmax : 1.9e-9 / Qh / seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 1,
}

def print_params(param_vals=param_vals):
    for i in param_vals:
        print(i, f' = {param_vals[i]:.2e}')

# functions
Xp = Bp / Qp
Xh = Bh / Qh

# vmax = muinfp* VmaxIp / Qp
Iuptakep = VmaxIp * (DIM / (DIM + KIp)) * exp(-omega*T) * Bp
Ouptakep = VmaxOp * (DOM / (DOM + KOp)) * exp(-omega*T) * Bp
Iuptakeh = VmaxIh * (DIM / (DIM + KIh)) * Bh
Ouptakeh = VmaxOh * (DOM / (DOM + KOh)) * Bh


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

dBpdt = Iuptakep + Ouptakep - deathp - exudationOp - exudationIp
dBhdt = Iuptakeh + Ouptakeh - deathh - exudationOh - exudationIh
dDOMdt = deathp * gammaDp + deathh * gammaDh + exudationOp +  exudationOh - Ouptakep - Ouptakeh
dRDOMdt = deathp * (1 - gammaDp) + deathh * (1 - gammaDh)
dDIMdt = exudationIp +  exudationIh - Iuptakep - Iuptakeh
dTdt = Treleasep - Tbreakdownh 

dDOMdt_ponly = deathp * gammaDp + exudationOp  - Ouptakep 
dRDOMdt_ponly = deathp * (1 - gammaDp) 
dDIMdt_ponly = exudationIp  - Iuptakep 
dTdt_ponly = Treleasep  

dDOMdt_honly = deathh * gammaDh +  exudationOh - Ouptakeh
dRDOMdt_honly = deathh * (1 - gammaDh)
dDIMdt_honly =  exudationIh -  Iuptakeh
dTdt_honly = - Tbreakdownh 


def print_equations():
    print('dBp/dt')
    display(dBpdt)
    print('dBh/dt')
    display(dBhdt)
    print('dDOM/dt')
    display(dDOMdt)
    print('dDIM/dt')
    display(dDIMdt)
    print('dT/dt')
    display(dTdt)

def get_main_data(param_vals=param_vals):
    sfunc_list = [dBpdt, dBhdt, dDOMdt, dRDOMdt, dDIMdt, dTdt]
    var_names = ['Bp', 'Bh',  'DOM', 'RDOM', 'DIM', 'T']
    init_vars = [1e9*Qp, 1e10*Qh, 20, 0, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bp, Bh, DOM, RDOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xp, Xh,
        Iuptakep, Ouptakep, Iuptakeh, Ouptakeh,
        deathp , deathh ,
        exudationOp, exudationIp, exudationOh, exudationIh, 
        Treleasep, Tbreakdownh,
    ]
    interm_names = [
        'Xp', 'Xh',
        'Iuptakep', 'Ouptakep', 'Iuptakeh', 'Ouptakeh',
        'deathp' , 'deathh' ,
        'exudationOp', 'exudationIp', 'exudationOh', 'exudationIh', 
        'Treleasep', 'Tbreakdownh',    
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bp, Bh, DOM, RDOM, DIM, T), interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_data(param_vals=param_vals):
    sfunc_list = [dBpdt, dDOMdt_ponly, dRDOMdt_ponly, dDIMdt_ponly, dTdt_ponly]
    var_names = ['Bp', 'DOM', 'RDOM', 'DIM', 'T']
    init_vars = [1e9*Qp, 20, 0, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bp, DOM, RDOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xp , 
        Iuptakep, Ouptakep, 
        deathp , 
        exudationOp, exudationIp, 
        Treleasep, 
    ]
    interm_names = [
        'Xp' , 
        'Iuptakep', 'Ouptakep', 
        'deathp' , 
        'exudationOp', 'exudationIp', 
        'Treleasep',  
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bp, DOM, RDOM, DIM, T), interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals=param_vals):
    sfunc_list = [dBhdt, dDOMdt_honly,  dRDOMdt_honly, dDIMdt_honly, dTdt_honly]
    var_names = ['Bh', 'DOM', 'RDOM', 'DIM', 'T']
    init_vars = [1e10*Qh, 20, 0, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bh, DOM, RDOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Xh , 
        Iuptakeh, Ouptakeh, 
        deathh , 
        exudationOh, exudationIh, 
        Tbreakdownh, 
    ]
    interm_names = [
        'Xh' , 
        'Iuptakeh', 'Ouptakeh', 
        'deathh' , 
        'exudationOh', 'exudationIh', 
        'Tbreakdownh',  
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bh, DOM, RDOM, DIM, T), interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func


def print_dydt0(calc_dydt, var_names, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    for i,j, k in zip(var_names, dydt0, init_vars):
        print(f'd{i}/dt = {j:.2e}, init {i} = {k:.2e}, newval = {k+j:.2e}')


def print_intermediate0(intermediate_func, interm_names, init_vars):
    for i,j in zip(interm_names, intermediate_func(*init_vars)):
        print(f'{i:<4} = {j:.2e}')


def biomass_diff0(calc_dydt, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    _Bp, _Bh, _DOM, _RDOM, _DIM, _T = dydt0
    print (f'dBp/dt + dBh/dt + dDOM/dt dRDOM/dt + dDIM/dt = {_Bp + _Bh + _DOM + _RDOM + _DIM}')

def biomass_diff0_ponly(calc_dydt, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    _Bp, _DOM, _RDOM, _DIM, _T = dydt0
    print (f'dBp/dt + dDOM/dt + dRDOM/dt + dDIM/dt = {_Bp + _DOM + _DIM + _RDOM}')

def biomass_diff0_honly(calc_dydt, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    _Bh, _DOM, _RDOM, _DIM, _T = dydt0
    print (f'dBh/dt + dDOM/dt + dRDOM/dt + dDIM/dt = {_Bh + _DOM + _DIM + _RDOM}')


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


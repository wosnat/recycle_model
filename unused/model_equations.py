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
Xp, Xh, Bp, Bh, DOM, DIM, T = symbols('X_p X_h B_p B_h DOM DIM T')

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

param_vals = {
    # 1/d
    Mh: 0.1/ seconds_in_day,
    Mp : 0.1/ seconds_in_day,
    # ratio
    gammaDp : 0.5,         # pro death release 
    gammaDh : 0.5,         # het death release
    # 1/d
    EOp : 0.2 / seconds_in_day,        # pro organic exudation 
    EIp : 0 / seconds_in_day,          # pro inorganic exudation
    EOh : 0 / seconds_in_day,          # het organic exudation 
    EIh : 0.2 / seconds_in_day,        # het inorganic exudation 
    # umol/cell
    Qmaxp : 1.5e-9, 
    Qminp : 7e-10, 
    Qmaxh : 1.5e-9 * pro_alt_vol_ratio,
    Qminh : 7e-10  * pro_alt_vol_ratio, 
    # TODO - change k for organic/inorganic
    # umol/l
    KOp : 0.17 * pro_vol**0.27, 
    KIp : 0.17 * pro_vol**0.27, 
    KOh : 0.17 * alt_vol**0.27, 
    KIh : 0.17 * alt_vol**0.27, 
    # umol/cell/d
    VmaxOp : 0, 
    VmaxIp : 1.9e-9 / seconds_in_day, 
    VmaxOh : 1.9e-9* pro_alt_vol_ratio**0.67 / seconds_in_day, 
    VmaxIh : 1.9e-9 / 10* pro_alt_vol_ratio**0.67 / seconds_in_day,
    # 1/d
    muinfp :  1 / seconds_in_day, 
    muinfh: 1.5 / seconds_in_day, 
    # umol/cell/d
    epsilon : 1e-10/ seconds_in_day, 
    VTmax : 1.9e-9/ seconds_in_day, 
    # umol/l
    KTh : 0.17 * alt_vol**0.27, 
    # 1/ umol/l
    omega : 1,
}

def print_params(param_vals=param_vals):
    for i in param_vals:
        print(i, f' = {param_vals[i]:.2e}')

# functions
# Qh = Max(Qminh, Bh / Xh)
# Qp = Max(Qminp, Bp / Xp)
Qh = Bh / Xh
Qp = Bp / Xp
Iuptakep = VmaxIp * (Qmaxp - Qp)/(Qmaxp - Qminp) * DIM / (DIM + KIp) * Xp
Ouptakep = VmaxOp * (Qmaxp - Qp)/(Qmaxp - Qminp) * DOM / (DOM + KOp) * Xp
Iuptakeh = VmaxIh * (Qmaxh - Qh)/(Qmaxh - Qminh) * DIM / (DIM + KIh) * Xh
Ouptakeh = VmaxOh * (Qmaxh - Qh)/(Qmaxh - Qminh) * DOM / (DOM + KOh) * Xh

mup  = muinfp* (1 - (Qminp / Qp)) * exp(-omega*T)
muh  = muinfh* (1 - (Qminh / Qh)) 
growthp = mup * Xp
growthh = muh * Xh


# if Q < Qmin : dX/dt = XX 
# newB/newX = Qmin
# newX = B / Qmin = B / Qmin 
# dX/dt = newX - X = B / Qmin - X = B /Qmin - B/ Q = B/(Qmin-Q)


Mstarp = Mp + Max(((Qminp - Qp)/Qminp), 0) / seconds_in_day
Mstarh = Mh + Max(((Qminh - Qh)/Qminh), 0) / seconds_in_day
deathp = Mstarp * Xp #* Xp
deathh = Mstarh * Xh #* Xh
exudationOp = EOp * Bp
exudationIp = EIp * Bp
exudationOh = EOh * Bh
exudationIh = EIh * Bh

Treleasep = epsilon * Xp
Tbreakdownh = VTmax * T / (T + KTh) * Xh

dBpdt = Iuptakep + Ouptakep - Qp * deathp - exudationOp - exudationIp
dBhdt = Iuptakeh + Ouptakeh - Qh * deathh - exudationOh - exudationIh
dXpdt = growthp - deathp
dXhdt = growthh - deathh
dDOMdt = Qp * deathp * gammaDp + Qh * deathh * gammaDh + exudationOp +  exudationOh - Ouptakep - Ouptakeh
dDIMdt = exudationIp +  exudationIh - Iuptakep - Iuptakeh
dTdt = Treleasep - Tbreakdownh 

dDOMdt_ponly = Qp * deathp * gammaDp + exudationOp  - Ouptakep 
dDIMdt_ponly = exudationIp  - Iuptakep 
dTdt_ponly = Treleasep  

dDOMdt_honly = Qh * deathh * gammaDh +  exudationOh - Ouptakeh
dDIMdt_honly =  exudationIh -  Iuptakeh
dTdt_honly = - Tbreakdownh 


def print_equations():
    print('dBp/dt')
    display(dBpdt)
    print('dBh/dt')
    display(dBhdt)
    print('dXp/dt')
    display(dXpdt)
    print('dXh/dt')
    display(dXhdt)
    print('dDOM/dt')
    display(dDOMdt)
    print('dDIM/dt')
    display(dDIMdt)
    print('dT/dt')
    display(dTdt)

def get_main_data(param_vals=param_vals):
    sfunc_list = [dBpdt, dBhdt, dXpdt, dXhdt, dDOMdt, dDIMdt, dTdt]
    var_names = ['Bp', 'Bh', 'Xp', 'Xh', 'DOM', 'DIM', 'T']
    init_vars = [(1e9+1e4)*param_vals[Qminp], (1e10+1e4)*param_vals[Qminh], 1e9, 1e10, 5, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bp, Bh, Xp, Xh, DOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Qh, Qp , 
        Iuptakep, Ouptakep, Iuptakeh, Ouptakeh,
        mup, muh, Mstarp, Mstarh,
        growthp, growthh,
        deathp , deathh ,
        exudationOp, exudationIp, exudationOh, exudationIh, 
        Treleasep, Tbreakdownh,
    ]
    interm_names = [
        'Qh', 'Qp' , 
        'Iuptakep', 'Ouptakep', 'Iuptakeh', 'Ouptakeh',
        'mup', 'muh', 'Mstarp', 'Mstarh',
        'growthp', 'growthh',
        'deathp' , 'deathh' ,
        'exudationOp', 'exudationIp', 'exudationOh', 'exudationIh', 
        'Treleasep', 'Tbreakdownh',    
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bp, Bh, Xp, Xh, DOM, DIM, T), interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_ponly_data(param_vals=param_vals):
    sfunc_list = [dBpdt, dXpdt, dDOMdt_ponly, dDIMdt_ponly, dTdt_ponly]
    var_names = ['Bp', 'Xp', 'DOM', 'DIM', 'T']
    init_vars = [(1e9+1e4)*param_vals[Qminp], 1e9, 5, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bp, Xp, DOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Qp , 
        Iuptakep, Ouptakep, 
        mup, Mstarp,
        growthp, 
        deathp , 
        exudationOp, exudationIp, 
        Treleasep, 
    ]
    interm_names = [
        'Qp' , 
        'Iuptakep', 'Ouptakep', 
        'mup', 'Mstarp',
        'growthp', 
        'deathp' , 
        'exudationOp', 'exudationIp', 
        'Treleasep',  
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bp, Xp, DOM, DIM, T), interm_funclist)

    return var_names, init_vars, calc_dydt, interm_names, intermediate_func

def get_honly_data(param_vals=param_vals):
    sfunc_list = [dBhdt, dXhdt, dDOMdt_honly, dDIMdt_honly, dTdt_honly]
    var_names = ['Bh', 'Xh', 'DOM', 'DIM', 'T']
    init_vars = [(1e10+1e4)*param_vals[Qminp], 1e10, 5, 100, 0] 
    subs_funclist = [sfunc.subs(param_vals) for sfunc in sfunc_list]
    final_func = lambdify((Bh, Xh, DOM, DIM, T), subs_funclist)
    calc_dydt = lambda t, y : final_func(*y)

    interm_sfunc_list = [
        Qh , 
        Iuptakeh, Ouptakeh, 
        muh, Mstarh,
        growthh, 
        deathh , 
        exudationOh, exudationIh, 
        Tbreakdownh, 
    ]
    interm_names = [
        'Qh' , 
        'Iuptakeh', 'Ouptakeh', 
        'muh', 'Mstarh',
        'growthh', 
        'deathh' , 
        'exudationOh', 'exudationIh', 
        'Tbreakdownh',  
    ]

    interm_funclist = [sfunc.subs(param_vals) for sfunc in interm_sfunc_list]
    intermediate_func = lambdify((Bh, Xh, DOM, DIM, T), interm_funclist)

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
    _Bp, _Bh, _Xp, _Xh, _DOM, _DIM, _T = dydt0
    print (f'dBp/dt + dBh/dt + dDOM/dt + dDIM/dt = {_Bp + _Bh + _DOM + _DIM}')

def biomass_diff0_ponly(calc_dydt, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    _Bp, _Xp, _DOM, _DIM, _T = dydt0
    print (f'dBp/dt + dDOM/dt + dDIM/dt = {_Bp + _DOM + _DIM}')

def biomass_diff0_honly(calc_dydt, init_vars):
    dydt0 = calc_dydt(0, init_vars)
    _Bh, _Xh, _DOM, _DIM, _T = dydt0
    print (f'dBh/dt + dDOM/dt + dDIM/dt = {_Bh + _DOM + _DIM}')


def run_solver(calc_dydt, init_vars, days=62, t_eval=None):
    tstart = 0
    tend = days*seconds_in_day
    if t_eval is None: 
        t_eval = np.arange(tstart, tend, 3600*4)
    sol = solve_ivp(
        fun=calc_dydt, y0=init_vars,
        t_span=[tstart, tend], t_eval=t_eval, max_step=100, first_step=1)
    print(f'solve_ivp(fun=calc_dydt, y0={init_vars},\n    t_span=[{tstart}, {tend}],\n    t_eval=t_eval)')
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


#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# %config InlineBackend.figure_formats = ['svg']


_new_black = '#373737'
plt_rc_params={
    'savefig.format' : 'svg', 
    'savefig.transparent' : True,
    #'font.family': 'sans-serif',
    #'font.sans-serif': ['Arial', 'DejaVu Sans'],
    # 'figure.dpi' : 600,
    'savefig.dpi' : 600,
    
    'legend.numpoints' : 1,
    'legend.markerscale' : 2.5,

    # 'figure.figsize': (15,8),
    
    'svg.fonttype': 'none',
    'text.usetex': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    # 'font.size': 9,
    # 'axes.labelsize': 9,
    # 'axes.titlesize': 9,
    # 'axes.labelpad': 2,
    # 'axes.linewidth': 0.5,
    # 'axes.titlepad': 4,
    # 'lines.linewidth': 0.5,
    # 'legend.fontsize': 9,
    # 'legend.title_fontsize': 9,
    # 'xtick.labelsize': 9,
    # 'ytick.labelsize': 9,
    # 'xtick.major.size': 2,
    # 'xtick.major.pad': 1,
    # 'xtick.major.width': 0.5,
    # 'ytick.major.size': 2,
    # 'ytick.major.pad': 1,
    # 'ytick.major.width': 0.5,
    # 'xtick.minor.size': 2,
    # 'xtick.minor.pad': 1,
    # 'xtick.minor.width': 0.5,
    # 'ytick.minor.size': 2,
    # 'ytick.minor.pad': 1,
    # 'ytick.minor.width': 0.5,

    # Avoid black unless necessary
    'text.color': _new_black,
    'patch.edgecolor': _new_black,
    'patch.force_edgecolor': False, # Seaborn turns on edgecolors for histograms by default and I don't like it
    'hatch.color': _new_black,
    'axes.edgecolor': _new_black,
    # 'axes.titlecolor': _new_black # should fallback to text.color
    'axes.labelcolor': _new_black,
    'xtick.color': _new_black,
    'ytick.color': _new_black

    # Default colormap - personal preference
    # 'image.cmap': 'inferno'
}
sns.set_style('white', rc=plt_rc_params)


gorder = [ 'Strong', 'Sustained', 'Inhibited', 'Weak', 'Neutral',   'Other']

gpalette = ['#882255', '#CC6677', '#332288', 
             '#44AA99','#88CCEE',
            '#D0CFCA',  ]
gorder1 = gorder[:-2]
gpalette = gpalette[:-2]
# sns.color_palette(gpalette)



axenic = ['Axenic']

group1 = ['A. macleodii 1A3', 'Pseudoalteromonas haloplanktis',]

group2 = ['Sulfitobacter pseudonitzschiae','Ruegeria pomeroyi', ]
group3 = [ 'Marinovum 5F3','Roseovarius 5C3']
group4 = [ 
       'Marinobacter adhaerens HP15',
       'Phaeobacter gallaeciensis', ]

horder =  group1 + group2 + group3 + group4 + ['Axenic']

# order based on tree order

horder_tree = ['Pseudoalteromonas haloplanktis',
 'A. macleodii 1A3',
 'Marinobacter adhaerens HP15',
 'Phaeobacter gallaeciensis',
 'Ruegeria pomeroyi',
 'Marinovum 5F3',
 'Roseovarius 5C3',
 'Sulfitobacter pseudonitzschiae']

gorder_tree = [
'Strong',
 'Strong',
 'Inhibited',
 'Inhibited',
 'Sustained',
 'Weak',
 'Weak',
 'Sustained']
gorder_full = [
'Strong',
 'Strong',
 'Sustained',
 'Sustained',
 'Weak',
 'Weak',
 'Inhibited',
 'Inhibited',
]

gpalette_dict = dict(zip(gorder, gpalette))
hpalette_tree = [gpalette_dict[i] for i in gorder_tree]
hpalette_g = [gpalette_dict[i] for i in gorder_full]
#sns.color_palette(hpalette_tree)
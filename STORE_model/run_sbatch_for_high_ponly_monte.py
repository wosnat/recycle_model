#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import re

def print_sbatch_file(model, vpro, template, outdpath):
    out_fpath = os.path.join(outdpath, f'run_monte_ponly_high_{model}_{vpro}.sbatch')
    with open(template) as vfh:
        with open(out_fpath, 'w') as ofh:
            for l in vfh:
                l =  re.sub(r'<MODEL>', model, l)
                l =  re.sub(r'<VPRO>', vpro, l)
                ofh.write(l)
    



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='post process results')
    # parser.add_argument("--model", help="model to run", choices=['MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME'], required=True)
    # parser.add_argument("--vprocsv", help="file with a list of VPROs to run",  required=True)
    parser.add_argument("--template", help="template file", default='run_files/hive/run_monte_ponly_high_template.sbatch')
    parser.add_argument("--outdpath", help="output dir", default='run_files/hive/run_monte_ponly_high/')
    

    args = parser.parse_args()
    outdpath = args.outdpath
    if outdpath != '':
        os.makedirs(outdpath, exist_ok=True)

    #for model in ['MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME']:
    for model in [ 'ROS', ]:
        #vprocsv = f'VPRO/high_growing_vpro_{model}.csv'
        vprocsv = f'VPRO/high_growing_vpro_round2_{model}.csv'
        with open(vprocsv) as vfh:
            for line in vfh:
                vpro = line.strip()
                print_sbatch_file(model, vpro, args.template, args.outdpath)


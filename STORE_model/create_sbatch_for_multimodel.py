#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import re

def print_sbatch_file(model, template, outdpath, vpro=None):
    templatebasename = os.path.basename(template)
    templatebasename = templatebasename.replace('.template','')
    out_fpath = os.path.join(outdpath, f'{templatebasename}_{model}.sbatch')
    if vpro is not None:
        vprobasename = os.path.basename(vpro)
        vprobasename = vprobasename.replace('.json','')
        out_fpath = os.path.join(outdpath, f'{templatebasename}_{model}_{vprobasename}.sbatch')

    with open(template) as vfh:
        with open(out_fpath, 'w') as ofh:
            for l in vfh:
                l =  re.sub(r'<MODEL>', model, l)
                if vpro is not None:
                    l =  re.sub(r'<VPRO>', vpro, l)
                ofh.write(l)
    



if __name__ == '__main__':
    import argparse

    multimodel_list = [
            'MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME',
            'ROS-MIXOTROPH-OVERFLOW-EXOENZYME', 'ROS-MIXOTROPH', 'EXOENZYME-ROS',
            'EXOENZYME-MIXOTROPH', 'OVERFLOW-MIXOTROPH', 'OVERFLOW-EXOENZYME', 'OVERFLOW-ROS',
            'EXOENZYME-ROS-MIXOTROPH', 'OVERFLOW-EXOENZYME-MIXOTROPH', 
            'OVERFLOW-ROS-MIXOTROPH', 'OVERFLOW-ROS-EXOENZYME',
        ]

    parser = argparse.ArgumentParser(description='create batch file for the multi models')
    #parser.add_argument("--model", help="model to run", choices=[
            #'MIN', 'MIXOTROPH', 'OVERFLOW', 'ROS', 'EXOENZYME',
            #'ROS-MIXOTROPH-OVERFLOW-EXOENZYME', 'ROS-MIXOTROPH', 'EXOENZYME-ROS',
            #'EXOENZYME-MIXOTROPH', 'OVERFLOW-MIXOTROPH', 'OVERFLOW-EXOENZYME', 'OVERFLOW-ROS',
            #'EXOENZYME-ROS-MIXOTROPH', 'OVERFLOW-EXOENZYME-MIXOTROPH', 
            #'OVERFLOW-ROS-MIXOTROPH', 'OVERFLOW-ROS-EXOENZYME',
        #], required=True)
    parser.add_argument("--template", help="template file", required=True)
    parser.add_argument("--outdpath", help="output dir", required=True)
    parser.add_argument("--vprodpath", help="vpro dir", default=None)
    

    args = parser.parse_args()
    outdpath = args.outdpath
    if outdpath != '':
        os.makedirs(outdpath, exist_ok=True)

    if args.vprodpath is not None:
        for model in multimodel_list:
            m_vpro_dpath = os.path.join(args.vprodpath, model)
            vprolst = os.listdir(m_vpro_dpath)
            for v in vprolst:
                vpro = os.path.join(m_vpro_dpath, v)
                print_sbatch_file(model, args.template, args.outdpath, vpro)
    else:
        for model in multimodel_list:
            print_sbatch_file(model, args.template, args.outdpath)


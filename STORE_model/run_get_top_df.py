#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import math
import os



if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='create df with only the top runs')

    parser.add_argument("--outdpath", help="output dir", default='.')
    parser.add_argument("--topids", help="csv file with a list of top runs", required=True)
    parser.add_argument("--infpath", help="input dir", default='.')
    
    
    args = parser.parse_args()
    out_dpath = args.outdpath
    if out_dpath != '':
        os.makedirs(out_dpath, exist_ok=True)

#######################

    fpath = args.infpath
    # results/final/het/monte_het_add_100per_vpro_ROS_clean_df.csv.gz
    fname = os.path.basename(fpath)
    session_id = fname.replace('_df.csv.gz', '')
    
    top_ids_df = pd.read_csv(args.topids)
    df = pd.read_csv(fpath)
    df = df.reset_index(drop=True)
    sample_df = df.loc[df.run_id.isin(top_ids_df.run_id)]

    sample_df.to_csv(os.path.join(out_dpath, f'top_runs_{session_id}_df.csv.gz'), index=False)

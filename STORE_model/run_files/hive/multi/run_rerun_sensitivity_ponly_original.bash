#!/bin/bash


RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR_base=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity_ponly
phase=rerun

fpath=$RDIR/results/final/ponly/vpros_init_sum.csv


    baseid=rerun_ponly_Original
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath"

#!/bin/bash


RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR_base=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity
phase=rerun


dpath=$RDIR/results/final/het

    # "monte_add_het_clean"
    # "monte_het_add_100per_vpro_EXOENZYME"
    # "monte_het_add_100per_vpro_OVERFLOW"
    # "monte_het_add_100per_vpro_ROS"
    # "monte_het_multi"
    # "monte_het_round2_100per_vpro_ROS"
    # "monte_ROS_round2_het"
fnames_sum=(
    "monte_het_clean"
    "monte_het_extend_100per_vpro_EXOENZYME"
    "monte_het_extend_100per_vpro_MIXOTROPH"
    "monte_het_extend_100per_vpro_OVERFLOW"
    "monte_het_extend_100per_vpro_ROS"
    "monte_het_extend_EXOENZYME"
    "monte_het_extend_MIXOTROPH"
    "monte_het_extend_OVERFLOW"
    "monte_het_extend_ROS"
)



for fname in "${fnames_sum[@]}"
do
    echo "Processing $fname"
    fpath="${dpath}/${fname}_clean_sum.csv.gz"
    baseid=rerun_het_DOC2500
    ODIR=${ODIR_base}/${baseid}/${fname}
    run_id=${baseid}
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism all --rerun_csv $fpath --override_init DOC,2500"
done



#!/bin/bash


RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR_base=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity_ponly
phase=rerun

fpath=$RDIR/results/final/ponly/vpros_init_sum.csv


    baseid=rerun_ponly_DIN800
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DIN,800"
    baseid=rerun_ponly_DIN800DOC160
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DIN,800 DOC,160"
    baseid=rerun_ponly_DIN800DOC2500
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DIN,800 DOC,2500"
    baseid=rerun_ponly_DIN800DOC2750
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DIN,800 DOC,2750"
    baseid=rerun_ponly_DOC160
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DOC,160"
    baseid=rerun_ponly_DOC2500
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DOC,2500"
    baseid=rerun_ponly_DOC2750
    ODIR=${ODIR_base}/${baseid}
    run_id=${baseid}
    echo $baseid
    sbatch --partition=hive7d,hiveunlim --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv $RDIR/reference_final.xlsx  --outdpath ${ODIR}/out --model MIN --run_id $run_id --which_organism ponly --rerun_csv $fpath --override_init DOC,2750"

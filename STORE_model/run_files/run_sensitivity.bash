
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/sensitivity_vmax2
run_id=sensitivity_vmax2

mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
    for o in HET
    do
       for i in  1 
       do
	       echo $m $o $i
	       $RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx  --outdpath ${ODIR}/out --run_id ${run_id} --model $m --organism_to_tune $o --param_sensitivity $i --number_of_runs 30 --which_organism all > ${ODIR}/log/${run_id}_${m}_ponly_${o}_${i}.log 2>&1 &



   done
   done
   done




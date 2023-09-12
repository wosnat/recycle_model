
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/sensitivity
run_id=sensitivity

mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
    for o in PRO
    do
       for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 
       do
	       echo $m $o $i
	       $RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx  --outdpath ${ODIR}/out --run_id ${run_id} --model $m --organism_to_tune $o --param_sensitivity $i --number_of_runs 20 --which_organism ponly > ${ODIR}/log/${run_id}_${m}_ponly_${o}_${i}.log 2>&1 &



   done
   sleep 200
   done
   done




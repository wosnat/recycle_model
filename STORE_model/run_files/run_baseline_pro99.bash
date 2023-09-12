
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/baseline
run_id=baseline

mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
    for w in ponly all
    do
	       echo $m $w 
	       run_id=baseline
	       $RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_pro99_axenic.xlsx  --pro99_mode --outdpath ${ODIR}/out --run_id ${run_id} --model $m --which_organism $w > ${ODIR}/log/${run_id}_${m}_${w}_p99.log 2>&1 &



   done
   done




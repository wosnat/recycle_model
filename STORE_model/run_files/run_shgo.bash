baseid=shgo
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
               echo $m 
               run_id=${baseid}
	       $RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx --ref_pro99_csv ${RDIR}/reference_pro99_axenic.xlsx  --outdpath ${ODIR}/out --run_id ${run_id} --model $m --which_organism $w --sobol 10 > ${ODIR}/log/${run_id}_${m}_${w}.log 2>&1 &




   done


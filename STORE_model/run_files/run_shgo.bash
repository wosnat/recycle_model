baseid=shgo
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
               echo $m 
               run_id=${baseid}_${m}_1
	       $RDIR/run_shgo_ponly.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx --ref_pro99_csv ${RDIR}/reference_pro99_axenic.xlsx   --out_dpath ${ODIR}/out --run_id ${run_id} --model $m --number_of_runs 10 > ${ODIR}/log/${run_id}.log 2>&1 &




   done


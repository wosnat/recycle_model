m=MIN
baseid=lsq_ponly_${m}
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/${baseid}
w=ponly


#MIN MIN OVERFLOW ROS MIN
#

mkdir -p $ODIR/log

echo $m 
for j in $(ls VPRO/X0/*${m}*)
do 
run_id=${baseid}_$(basename ${j%.*})
echo $m $j $run_id
$RDIR/run_least_squares_ponly.py --ref_csv $RDIR/reference_10cc_axenic.xlsx --ref_pro99_csv $RDIR/reference_pro99_axenic.xlsx --out_dpath ${ODIR}/out --run_id ${run_id} --model $m --json $j > $ODIR/log/${run_id}.log 2>&1 &



#~/work/RECYCLE_MODEL/recycle_model/model_equations_separate_NC_sep_vmax.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc_axenic.xlsx --outdpath ~/work/RECYCLE_MODEL/results/vpro/out/ --run_id vpro_ponly_${run_id} --model ${m}  --which_organism ponly --json ${j} > ~/work/RECYCLE_MODEL/results/vpro/log/${run_id}.log 2>&1 &


   done








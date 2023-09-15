baseid=vpro_ponly
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}

w=ponly
mkdir -p $ODIR/log

for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do 
for j in $(ls $PWD/VPRO/FINAL/*${m}*.json)
do 
vpro=${j##*/}
vpro=${vpro%.json}
echo $m $vpro
run_id=${baseid}_${vpro}
sbatch --job-name=${run_id} --partition=hive1d,hive7d,hiveunlim --output=${ODIR}/log/out_%A_%j.out  --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx  --outdpath ${ODIR}/out --run_id ${run_id} --model $m --which_organism $w --json $j"
sbatch --job-name=${run_id}_p99 --partition=hive1d,hive7d,hiveunlim --output=${ODIR}/log/out_%A_%j.out  --wrap "$RDIR/model_equations_separate_NC_store_numba.py --ref_csv ${RDIR}/reference_pro99_axenic.xlsx --pro99_mode  --outdpath ${ODIR}/out --run_id ${run_id} --model $m --which_organism $w --json $j"

done; 
done











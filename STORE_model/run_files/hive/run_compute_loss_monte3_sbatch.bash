#!/bin/bash

#SBATCH --job-name=lsqp_MIN
#SBATCH --partition=hive1d,hive7d,hiveunlim
#SBATCH --array=11-50
#SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID



#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands

echo "branch${i}";






#MIN MIXOTROPH OVERFLOW ROS EXOENZYME
#
i=0
baseid=loss_monte3
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log
mkdir -p $ODIR/out
for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
for loss in linear huber softl1 cauchy arctan
do
for f_scale in 0.01 0.1 1 10 100
do
run_id=${baseid}_${m}_${loss}_${f_scale}
echo $run_id
sbatch --job-name=${run_id} --partition=hive1d,hive7d,hiveunlim --output=${ODIR}/log/${run_id}_%A_%j.out  --wrap "${RDIR}/run_compute_mce.py --ref_csv ${RDIR}/reference_10cc.xlsx --outdpath ${ODIR}/out --run_id ${run_id} --incsv $RDIR/results/monte3_${m}_df.csv.gz --which all --loss ${loss} --f_scale ${f_scale}"


done
done
done





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
baseid=lsq_ponly
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log
for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
do
echo $m 
run_id=${baseid}_${m}_${i}
sbatch --job-name=lsqp_${m} --partition=hive1d,hive7d,hiveunlim --output=${ODIR}/log/out_%A_%j.out  --wrap "${RDIR}/run_least_squares_ponly.py --ref_csv ${RDIR}/reference_10cc_axenic.xlsx --ref_pro99_csv ${RDIR}/reference_pro99_axenic.xlsx --out_dpath ${ODIR}/out --run_id ${run_id} --model $m"


done




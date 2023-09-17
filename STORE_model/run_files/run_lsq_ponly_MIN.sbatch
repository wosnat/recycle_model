#!/bin/bash

#SBATCH --job-name=lsqp_MIN
#SBATCH --partition=hive1d,hive7d,hiveunlim
#SBATCH --array=11-50
#SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID

i=${SLURM_ARRAY_TASK_ID}


#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands

echo "branch${i}";






#MIN MIN OVERFLOW ROS MIN
#
m=MIN
baseid=lsq_ponly
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log

echo $m 
run_id=${baseid}_${m}_${i}
$RDIR/run_least_squares_ponly.py --ref_csv $RDIR/reference_10cc_axenic.xlsx --ref_pro99_csv $RDIR/reference_pro99_axenic.xlsx --out_dpath ${ODIR}/out --run_id ${run_id} --model $m --json VPRO/X0/pro_x0_${m}_${i}.json







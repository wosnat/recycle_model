#!/bin/bash

#SBATCH --job-name=monte_MIN
#SBATCH --partition=hive1d,hive7d,hiveunlim
#SBATCH --array=1-300
# #SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID

i=${SLURM_ARRAY_TASK_ID}


#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands







#MIN MIN OVERFLOW ROS MIN
#
m=MIN
baseid=monte_extend_${m}
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/multi/het_extend/${baseid}
w=ponly

mkdir -p $ODIR/log

echo $m 
run_id=${baseid}_${i}
$RDIR/model_equations_separate_NC_store_numba.py --monte --ref_csv $RDIR/reference_final.xlsx  --jsondpath $RDIR/VPRO/multi/run2/${m}/ --outdpath ${ODIR}/out --model ${m} --run_id $run_id --which_organism all --organism_to_tune HET --number_of_runs 100




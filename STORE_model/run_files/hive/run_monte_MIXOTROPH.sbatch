#!/bin/bash

#SBATCH --job-name=monte_MIXOTROPH
#SBATCH --partition=hive1d,hive7d,hiveunlim
#SBATCH --array=1-100
#SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID

i=${SLURM_ARRAY_TASK_ID}


#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands

echo "branch${i}";






#MIXOTROPH MIXOTROPH OVERFLOW ROS MIXOTROPH
#
m=MIXOTROPH
baseid=monte
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
w=ponly

mkdir -p $ODIR/log

echo $m 
run_id=${baseid}_${m}_${i}
$RDIR/model_equations_separate_NC_store_numba.py --monte --ref_csv $RDIR/reference_10cc.xlsx  --jsondpath $RDIR/VPRO/LSQ/${m}/ --outdpath ${ODIR}/out --model ${m} --run_id $run_id --which_organism all --organism_to_tune HET --number_of_runs 100






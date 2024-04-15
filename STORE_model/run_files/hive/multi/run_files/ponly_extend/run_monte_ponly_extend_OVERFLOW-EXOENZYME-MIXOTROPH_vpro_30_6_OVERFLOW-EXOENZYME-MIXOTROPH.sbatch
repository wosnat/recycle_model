#!/bin/bash

#SBATCH --job-name=monte_OVERFLOW-EXOENZYME-MIXOTROPH
#SBATCH --partition=hive1d,hive7d,hiveunlim
# #SBATCH --array=1-100
#  #SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID

i=${SLURM_ARRAY_TASK_ID}


#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands

echo "branch${i}";






#EXOENZYME EXOENZYME OVERFLOW ROS EXOENZYME
#
m=OVERFLOW-EXOENZYME-MIXOTROPH
vpro=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/VPRO/multi/run1/OVERFLOW-EXOENZYME-MIXOTROPH/vpro_30_6_OVERFLOW-EXOENZYME-MIXOTROPH.json
vproname=$(basename $vpro .json)
baseid=monte_add_ponly_${m}_${vproname}
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/multi/ponly_extend/${baseid}
w=ponly

mkdir -p $ODIR/log

echo $m 
run_id=${baseid}_${i}

$RDIR/model_equations_separate_NC_store_numba.py --monte --ref_csv ${RDIR}/reference_10cc_axenic.xlsx --ref_pro99_csv ${RDIR}/reference_pro99_axenic.xlsx --outdpath ${ODIR}/out --run_id ${run_id} --model $m --which_organism ponly --organism_to_tune PRO --number_of_runs 100 --json ${vpro} --monte_max_params 5






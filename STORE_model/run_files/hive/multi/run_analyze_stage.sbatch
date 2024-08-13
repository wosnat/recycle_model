#!/bin/bash

#SBATCH --job-name=stage
#SBATCH --partition=hive1d,hive7d,hiveunlim
#SBATCH --array=0-39
# #SBATCH --output=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/out_%A_%a_%j.out

## To make things simple: ${i} == $SLURM_ARRAY_TASK_ID

i=${SLURM_ARRAY_TASK_ID}


#load the modules environment

#. /etc/profile.d/modules.sh
#module load Miniconda3

#Job commands



RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=${RDIR}/results/final/het/stages

filelist=( $(ls  ${RDIR}/results/final/het/*clean_df*) )
FPATH=${filelist[i]}

mkdir -p $ODIR

echo $RDIR/run_het_stage_df.py --outdpath $ODIR --infpath $FPATH
$RDIR/run_het_stage_df.py --outdpath $ODIR --infpath $FPATH





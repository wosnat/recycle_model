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
ODIR=${RDIR}/results/final/het/biogeo
dpath=${RDIR}/results/final/het

filelist=( $(ls  ${dpath}/*clean_df*) )

FPATH=${filelist[i]}
fname="${FPATH##*/}"
runid="${fname%_clean_df.csv.gz}"

mkdir -p $ODIR

echo $RDIR/run_het_biogeo_create_files.py --outdpath $ODIR --indpath $dpath --run_id $runid
$RDIR/run_het_biogeo_create_files.py --outdpath $ODIR --indpath $dpath --run_id $runid
echo done


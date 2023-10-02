#!/usr/bin/env bash






# run by: 
# cd slurm
# RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model

# for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done; done



# this script is meant to be run using gnu parallel.
# run using:
#  m=MIN; ls VPRO/X0/*${m}* |  parallel run_files/run_lsq_ponly_parallel.bash $m {}
# for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME
set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#  set -x  # Similar to verbose mode (-v), but expands commands
set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.



m=$1
j=$2

baseid=lsq_ponly_${m}
RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}

#RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
#ODIR=~/work/RECYCLE_MODEL/results/${baseid}
w=ponly


#MIN MIN OVERFLOW ROS MIN
#

mkdir -p $ODIR/log

echo $m 
#for j in $(ls VPRO/X0/*${m}*)
#do 
run_id=${baseid}_$(basename ${j%.*})
echo $m $j $run_id
$RDIR/run_least_squares_ponly.py --ref_csv $RDIR/reference_10cc_axenic.xlsx --ref_pro99_csv $RDIR/reference_pro99_axenic.xlsx --out_dpath ${ODIR}/out --run_id ${run_id} --model $m --json $j --logerror > $ODIR/log/${run_id}.log 2>&1 



#~/work/RECYCLE_MODEL/recycle_model/model_equations_separate_NC_sep_vmax.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc_axenic.xlsx --outdpath ~/work/RECYCLE_MODEL/results/vpro/out/ --run_id vpro_ponly_${run_id} --model ${m}  --which_organism ponly --json ${j} > ~/work/RECYCLE_MODEL/results/vpro/log/${run_id}.log 2>&1 &


#   done








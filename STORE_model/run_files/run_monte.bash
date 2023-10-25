#!/usr/bin/env bash


# this script is meant to be run using gnu parallel.
# run using:
#   for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do parallel run_files/run_monte.bash $m ::: {1..100}; done

set -e  # Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
set -u  # Attempt to use undefined variable outputs error message, and forces an exit
#  set -x  # Similar to verbose mode (-v), but expands commands
set -o pipefail  # Causes a pipeline to return the exit status of the last command in the pipe that returned a non-zero return value.



m=$1
i=$2

baseid=monte4
RDIR=/home/oweissberg/work/RECYCLE_MODEL/recycle_model/STORE_model
ODIR=~/work/RECYCLE_MODEL/results/${baseid}


#MIN MIN OVERFLOW ROS MIN
#

mkdir -p $ODIR/log

run_id=${baseid}_${m}_${i}
echo $m $i $run_id
$RDIR/model_equations_separate_NC_store_numba.py --monte --ref_csv $RDIR/reference_10cc.xlsx  --jsondpath $RDIR/VPRO/X0/${m}/ --outdpath ${ODIR}/out --model ${m} --run_id $run_id --which_organism all --organism_to_tune HET --number_of_runs 100









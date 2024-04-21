#!/bin/bash


RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
ODIR_base=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/multi/het_100per_vpro
phase=extend

# EXOENZYME
# EXOENZYME-MIXOTROPH
# EXOENZYME-ROS
# EXOENZYME-ROS-MIXOTROPH
# MIN
# MIXOTROPH
# OVERFLOW
# OVERFLOW-EXOENZYME
# OVERFLOW-EXOENZYME-MIXOTROPH
# OVERFLOW-MIXOTROPH
# OVERFLOW-ROS
# OVERFLOW-ROS-EXOENZYME
# OVERFLOW-ROS-MIXOTROPH
# ROS
# ROS-MIXOTROPH
# ROS-MIXOTROPH-OVERFLOW-EXOENZYME

# only those starting with EXO
prefix=${RDIR}/VPRO/multi/run2/
for d in ${prefix}EXO*
do 
	m=${d/$prefix/}
	echo $m
	for f in $d/*
	do
		vpro=${f/$d\//}
		vpro=${vpro/.json/}
		echo $d $m $vpro
		baseid=monte_rerun_het_${phase}
		ODIR=${ODIR_base}2/${phase}/${m}
		run_id=${baseid}_${vpro}
		sbatch --partition=hive1d,hive7d,hiveunlim,preempt1d,preempt7d,preempt31d --wrap "$RDIR/model_equations_separate_NC_store_numba.py --monte --ref_csv $RDIR/reference_final.xlsx  --json $f --outdpath ${ODIR}/out --model ${m} --run_id $run_id --which_organism all --organism_to_tune HET --number_of_runs 100"
	done
done



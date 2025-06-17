
RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do echo $i; sbatch --wrap "$PWD/results_postprocess.py --dpath ~/RECYCLE_MODEL/results/rerun_sensitivity/$i/*/out --run_id $i --outdpath results/rerun_sensitivity/"; done

RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do echo $i; sbatch --wrap "$PWD/run_het_clean_classify_veratiles.py --indpath results/rerun_sensitivity/ --run_id $i --outdpath results/rerun_sensitivity/clean"; done


# init only

ls -d ~/RECYCLE_MODEL/results/rerun_sensitivity/*/monte_het_clean/out

RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do echo $i; sbatch --wrap "$PWD/results_postprocess.py --dpath ~/RECYCLE_MODEL/results/rerun_sensitivity/$i/monte_het_clean/out --run_id ${i}_init --outdpath results/rerun_sensitivity_init/"; done

RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do i=${i}_init; echo $i; sbatch --wrap "$PWD/run_het_clean_classify_veratiles.py --indpath results/rerun_sensitivity_init/ --run_id $i --outdpath results/rerun_sensitivity_init/clean"; done

cd slurm
sbatch ../run_files/hive/multi/run_create_biogeo_rerun_sensitivity.sbatch


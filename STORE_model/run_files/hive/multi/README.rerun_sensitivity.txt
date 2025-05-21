
RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do echo $i; sbatch --wrap "$PWD/results_postprocess.py --dpath ~/RECYCLE_MODEL/results/rerun_sensitivity/$i/*/out --run_id $i --outdpath results/rerun_sensitivity/"; done

RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/rerun_sensitivity/; for i in $(ls $RDIR); do echo $i; sbatch --wrap "$PWD/run_het_clean_classify_veratiles.py --indpath results/rerun_sensitivity/ --run_id $i --outdpath results/rerun_sensitivity/clean"; done


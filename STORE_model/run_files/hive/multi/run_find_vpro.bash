RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model

#--outdpath OUTDPATH   output dir
#--vprooutdpath VPROOUTDPATH
#--prefixvpro PREFIXVPRO
#--indpath INDPATH     input dir
#--run_id RUN_ID       run id

dpath=$RDIR/results/multi/ponly

#mkdir -p ./results/multi/ponly;dp=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/multi/monte_ponly_;for i in ~/RECYCLE_MODEL/results/multi/monte_ponly_*; do m=${i/${dp}/}; sbatch --wrap "./results_postprocess.py --dpath $i/out/ --outdpath ./results/multi/ponly/ --run_id monte_ponly_${m}" ; done

dp=$dpath/monte_ponly_
for i in $dpath/monte_ponly_*_df.csv.gz
do 
	m=${i/${dp}/}; 
	m=${m/_df.csv.gz/}; 

	sbatch --wrap "../run_ponly_clean_find_vpros.py --outdpath $dpath --indpath $dpath --vprooutdpath $RDIR/VPRO/multi/run1/ --prefixvpro monte_ponly_${m} --run_id monte_ponly_${m}"

done

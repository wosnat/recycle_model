
RDIR=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model

#--outdpath OUTDPATH   output dir
#--vprooutdpath VPROOUTDPATH
#--prefixvpro PREFIXVPRO
#--indpath INDPATH     input dir
#--run_id RUN_ID       run id

dpath=$RDIR/results/multi/extend/het
runprefix=monte_het_extend_
#mkdir -p ./results/multi/ponly;dp=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/multi/monte_ponly_;for i in ~/RECYCLE_MODEL/results/multi/monte_ponly_*; do m=${i/${dp}/}; sbatch --wrap "./results_postprocess.py --dpath $i/out/ --outdpath ./results/multi/ponly/ --run_id monte_ponly_${m}" ; done

dp=$dpath/$runprefix
#for i in $dp*_df.csv.gz
for i in ${dp}{MIN,MIXOTROPH}_df.csv.gz
do 
	m=${i/${dp}/}; 
	m=${m/_df.csv.gz/}; 

	sbatch --wrap "$RDIR/run_het_clean_classify_veratiles.py --indpath $dpath --run_id ${runprefix}${m} --outdpath $dpath/clean"
done

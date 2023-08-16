

for m in MIXO
do
	       echo $m 
mkdir -p ~/work/RECYCLE_MODEL/results/shgo_${m}/log
~/work/RECYCLE_MODEL/recycle_model/run_shgo_ponly.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc_axenic.xlsx --ref_pro99_csv ~/work/RECYCLE_MODEL/recycle_model/reference_pro99_axenic.xlsx --json_dpath ~/work/RECYCLE_MODEL/results/shgo_${m}/json --out_dpath ~/work/RECYCLE_MODEL/results/shgo_${m}/out/ --run_id shgo_${m} --timeout 120 --model ${m}  > ~/work/RECYCLE_MODEL/results/shgo_${m}/log/shgo_${m}.log 2>&1 &



   done




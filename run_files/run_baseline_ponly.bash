
mkdir -p ~/work/RECYCLE_MODEL/results/baseline/log

for m in MIN LEAK MIXO FULL
do
	       echo $m 
~/work/RECYCLE_MODEL/recycle_model/model_equations_separate_NC_sep_vmax.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc_axenic.xlsx --json_dpath ~/work/RECYCLE_MODEL/results/baseline/json --out_dpath ~/work/RECYCLE_MODEL/results/baseline/out/ --run_id baseline_ponly --timeout 300 --model ${m}  --which_organism ponly > ~/work/RECYCLE_MODEL/results/baseline/log/baseline_ponly_${m}.log 2>&1 &


   done




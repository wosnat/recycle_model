mkdir -p ~/work/RECYCLE_MODEL/results/sensitivity_pro99/log

for m in MIN LEAK MIXO FULL
do
    for o in PRO
    do
       for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
       do
	       echo $m $o $i
~/work/RECYCLE_MODEL/recycle_model/run_sensitivity.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/refence_pro99_axenic.xlsx --pro99_mode --json_dpath ~/work/RECYCLE_MODEL/results/sensitivity_pro99/json --out_dpath ~/work/RECYCLE_MODEL/results/sensitivity_pro99/out/ --run_id param_sensitivity_ponly --timeout 300 --model ${m}  --organism_to_tune ${o} --param_sensitivity ${i} --number_of_runs 20 --which_organism ponly > ~/work/RECYCLE_MODEL/results/sensitivity_pro99/log/sensitivity_ponly_${m}_${o}_${i}.log 2>&1 &


   done
   sleep 1000
   done
   done






for m in MIN, LEAK, MIXO, FULL
do
    for o in PRO, HET
    do
       for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
       do
~/work/RECYCLE_MODEL/recycle_model/run_sensitivity.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/prelim_bottle.csv --json_dpath ~/work/RECYCLE_MODEL/results/sensitivity/json --out_dpath ~/work/RECYCLE_MODEL/results/sensitivity/out/ --run_id param_sensitivity_${m}_${o} --timeout 300 --model ${m}  --organism_to_tune ${o} --param_sensitivity ${i} --number_of_runs 20 > ~/work/RECYCLE_MODEL/results/sensitivity/log/sensitivity_${m}_${o}_${i}.log 2>&1 &


   done
   sleep 3600
   done
   done



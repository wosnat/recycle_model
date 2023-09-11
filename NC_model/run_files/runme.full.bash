
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
do
echo " ~/work/RECYCLE_MODEL/recycle_model/run_sensitivity.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/prepare_data/refdf_10cc_${i}.csv.gz --json_dpath ~/work/RECYCLE_MODEL/results/10cc/json/de10_full_${i} --out_dpath ~/work/RECYCLE_MODEL/results/10cc/out/de10_full_${i} --run_id de_10_full_lab_${i} --timeout 300 --optimize --model full > de10_full_${i}.log 2>&1 & "
# {min,full,exu,mix}

   done


for m in LEAK 
do
       for i in  0.1 0.01 1 10 20 100
       #for i in 1
       do
	       echo $m $i
mkdir -p ~/work/RECYCLE_MODEL/results/least_square_${m}_${i}/log
~/work/RECYCLE_MODEL/recycle_model/run_least_squares_ponly.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc_axenic.xlsx --ref_pro99_csv ~/work/RECYCLE_MODEL/recycle_model/reference_pro99_axenic.xlsx --json_dpath ~/work/RECYCLE_MODEL/results/least_square_${m}_${i}/json --out_dpath ~/work/RECYCLE_MODEL/results/least_square_${m}_${i}/out/ --run_id least_square_${m}_${i} --timeout 120 --model ${m} --f_scale ${i}  > ~/work/RECYCLE_MODEL/results/least_square_${m}_${i}/log/least_square_${m}_${i}.log 2>&1 &



   done
   done





# done upto including 6
for m in LEAK 
do
	pre_m=MIN
       #for i in  0.1 0.01  10 20 100
       for i in 1
       do
	       for s in 2
	       do
for j in $(ls VHET_JSON/*${pre_m}*_${s}_1.json)
do 
run_id=$(basename ${j%.*})
short_json=$(basename ${j%.*})
vpro_id=${short_json/het_least_square_MIN_/}; 
vpro_id1=${vpro_id/%_?_1/};
vpro_id2=${vpro_id1/%_1?_1/}; 

run_id=het_least_square_${m}_${short_json}_${s}_${i}
echo $m $i $s $j $run_id

mkdir -p ~/work/RECYCLE_MODEL/results/lsq_${m}/${run_id}/log


~/work/RECYCLE_MODEL/recycle_model/run_least_squares_het.py --ref_csv  ~/work/RECYCLE_MODEL/recycle_model/reference_10cc.xlsx --ref_sample_id ${s} --vpro_json PRO_JSON/${vpro_id2}.json ${j}  --json_dpath ~/work/RECYCLE_MODEL/results/lsq_${m}/${run_id}/json --out_dpath ~/work/RECYCLE_MODEL/results/lsq_${m}/${run_id}/out --run_id ${run_id} --model $m --f_scale $i --timeout 200  > ~/work/RECYCLE_MODEL/results/lsq_${m}/${run_id}/log/${run_id}.log 2>&1 &


   done
   done
   done
   done




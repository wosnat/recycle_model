

for m in MIXO
do
       #for i in  0.1 0.01  10 20 100
       for j in $(ls PRO_JSON/*_${m}*.json) 
       do

 short_json=$(basename ${j%.*})
 run_id=het_monte_${m}_${short_json}
	       echo $m  $j $run_id

mkdir -p ~/work/RECYCLE_MODEL/results/${run_id}/log


~/work/RECYCLE_MODEL/recycle_model/run_monte.py --ref_csv  ~/work/RECYCLE_MODEL/recycle_model/reference_10cc.xlsx --vpro_json ${j} --json_dpath ~/work/RECYCLE_MODEL/results/${run_id}/json --out_dpath ~/work/RECYCLE_MODEL/results/${run_id}/out --run_id ${run_id} --model $m  --number_of_runs 1024 --timeout 200  > ~/work/RECYCLE_MODEL/results/${run_id}/log/${run_id}.log 2>&1 &


   done
   done




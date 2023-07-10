
mkdir -p ~/work/RECYCLE_MODEL/results/vhet/log

for m in MIN 
do
for j in $(ls VHET_JSON/*${m}*)
do 
run_id=$(basename ${j%.*})
vpro_id=${run_id/het_least_square_MIN_/}; 
vpro_id1=${vpro_id/%_?_1/};
vpro_id2=${vpro_id1/%_1?_1/}; 

echo $m $j $run_id $vpro_id2
~/work/RECYCLE_MODEL/recycle_model/model_equations_separate_NC_sep_vmax.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_10cc.xlsx --outdpath ~/work/RECYCLE_MODEL/results/vhet/out/ --run_id vhet_${run_id} --model ${m}  --which_organism all --json PRO_JSON/${vpro_id2}.json ${j} > ~/work/RECYCLE_MODEL/results/vhet/log/${run_id}.log 2>&1 &

   done
   done




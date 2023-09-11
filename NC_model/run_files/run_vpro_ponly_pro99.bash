
mkdir -p ~/work/RECYCLE_MODEL/results/vpro/log

for m in MIN LEAK MIXO
do
for j in $(ls PRO_JSON/*${m}*)
do 
run_id=$(basename ${j%.*})
echo $m $j $run_id
~/work/RECYCLE_MODEL/recycle_model/model_equations_separate_NC_sep_vmax.py --ref_csv ~/work/RECYCLE_MODEL/recycle_model/reference_pro99_axenic.xlsx --pro99_mode --outdpath ~/work/RECYCLE_MODEL/results/vpro/out/ --run_id vpro_ponly_pro99_${run_id} --model ${m}  --which_organism ponly --json ${j} > ~/work/RECYCLE_MODEL/results/vpro/log/${run_id}_pro99.log 2>&1 &


   done
   done




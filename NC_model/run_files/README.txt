# to run sensitivity:

tmux new -s sensitivity

cd ~/work/RECYCLE_MODEL/recycle_model
source run_files/run_sensitivity.bash

# post process:
./results_postprocess.py --dpath ~/work/RECYCLE_MODEL/results/sensitivity/out/ --outdpath results --run_id param_sensitivity_29062023



./results_postprocess.py --dpath ~/work/RECYCLE_MODEL/results/baseline/out/ --outdpath results --run_id baseline_29062023

